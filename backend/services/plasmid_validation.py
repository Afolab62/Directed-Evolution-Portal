from __future__ import annotations

"""
Part C — Plasmid validation (circular DNA → ORF → protein)

Goal:
Validate that a *circular* plasmid DNA sequence encodes a given wild-type (WT) protein.

Approach (tiered, robust):
A) Exact match in 6-frame translation (supports 'X' wildcard for ambiguous DNA).
B) Fuzzy same-length window identity (handles substitutions; no indels).
C) Smith–Waterman local alignment (handles indels; slower, guarded by size limits).

Design notes:
- Plasmids are circular, so we search on dna2 = dna + dna and map coordinates back into [0, L).
- Coordinate mapping is defined on the *original plasmid*:
    start_nt:       0-based index of first coding nucleotide on original plasmid.
    end_nt_exclusive: 0-based exclusive end on original plasmid.
    wraps_origin:   True if the called ORF crosses the origin (i.e., end wraps past L).
  When wraps_origin=True, end_nt_exclusive will be < start_nt (modulo representation).
- We return rich diagnostics on failure so users/markers can see *why* validation failed.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Literal, Optional, Tuple

from .sequence_tools import CODON_TABLE, translate_six_frames

Strand = Literal["+", "-"]
MatchType = Literal["exact", "fuzzy", "align", "none"]


@dataclass(frozen=True)
class GeneCall:
    """
    Primary result object for plasmid validation.

    identity / coverage:
      - For exact/fuzzy: identity=coverage=1.0 (exact) or identity>=threshold, coverage=1.0 (fuzzy).
      - For alignment: identity and coverage are computed from the aligned query (WT) positions.

    diagnostics:
      - On failure: dict containing params + top candidate summaries.
      - On success: None by default; may contain plausibility_warnings if enabled.
    """
    is_valid: bool
    strand: Optional[Strand]
    frame: Optional[int]  # 0,1,2 for frame offset; strand indicates direction
    start_nt: Optional[int]
    end_nt_exclusive: Optional[int]
    wraps_origin: bool
    match_type: MatchType
    identity: Optional[float] = None
    coverage: Optional[float] = None
    notes: str = ""
    diagnostics: Optional[Dict[str, Any]] = None


@dataclass
class _Candidate:
    strand: Strand
    frame: int
    frame_key: str
    aa_index: int  # start in AA coordinates within the translated frame of dna2 (or rc(dna2))
    nt_start2: int  # start in nt coords on dna2 (0..2L)
    nt_end2: int    # end exclusive in nt coords on dna2 (0..2L)
    wraps_origin: bool
    score: int
    match_type: MatchType
    identity: float
    coverage: float
    notes: str


# ----------------------------
# Small utilities
# ----------------------------

def _clean_dna(text: str) -> str:
    # normalise to uppercase and strip whitespace/newlines — FASTA files
    # often have line breaks mid-sequence that would break codon reading
    return text.strip().upper().replace("\n", "").replace("\r", "")


def _clean_protein(text: str) -> str:
    return text.strip().upper().replace("\n", "").replace("\r", "")


def _frame_key_to_strand_frame(k: str) -> Tuple[Strand, int]:
    # translate_six_frames returns keys like "+0","+1","+2","-0","-1","-2"
    strand: Strand = "+" if k.startswith("+") else "-"
    frame = int(k[1:])
    return strand, frame


def _map_plus_nt_coords(L: int, nt_start2: int, length_nt: int) -> Tuple[int, int, bool]:
    """
    Map a match found on dna2 (length 2L) back to original plasmid (length L).
    Returns (start_nt, end_nt_exclusive, wraps_origin).

    Convention:
      - start_nt is modulo L
      - wraps_origin is True if the ORF crosses the origin on the original circle
      - if wraps_origin=True, end_nt_exclusive will be modulo'd so end < start
    """
    start_nt = nt_start2 % L
    nt_end2 = nt_start2 + length_nt

    # Wraps the origin if it starts in the first copy but ends beyond L
    wraps_origin = (nt_start2 < L) and (nt_end2 > L)

    if not wraps_origin:
        end_nt_exclusive = start_nt + length_nt
        # If start was in second copy, end could exceed L; keep within [0, L] by modulo only if needed.
        if end_nt_exclusive > L:
            end_nt_exclusive = end_nt_exclusive % L
    else:
        end_nt_exclusive = (start_nt + length_nt) - L  # end < start
    return start_nt, end_nt_exclusive, wraps_origin


def _map_minus_nt_coords(L: int, dna2_len: int, nt_start_rev2: int, length_nt: int) -> Tuple[int, int, bool, int]:
    """
    For '-' strand frames, the match is found on the reverse-complement of dna2.
    Convert start position in the rc string back to a start coordinate on original dna2 (plus strand).

    If rc index is nt_start_rev2, then the corresponding region on dna2 starts at:
      start_nt2 = dna2_len - (nt_start_rev2 + length_nt)

    Returns:
      (start_nt, end_nt_exclusive, wraps_origin, start_nt2) where start_nt2 is on dna2.
    """
    start_nt2 = dna2_len - (nt_start_rev2 + length_nt)
    start_nt, end_nt_exclusive, wraps_origin = _map_plus_nt_coords(L, start_nt2, length_nt)
    return start_nt, end_nt_exclusive, wraps_origin, start_nt2


def _orf_context_score(plasmid_dna: str, strand: Strand, start_nt: int, end_nt_exclusive: int) -> int:
    """
    Lightweight tie-breaker to prefer biologically plausible ORFs.
    This is *not* used to reject; only to choose between multiple equally-good matches.

    Heuristic:
      +5 if start codon at start is ATG (for '+' only)
      +5 if a canonical stop codon immediately follows end (for '+' only)
    For '-' we skip the codon checks (still tie-break stable via 0).

    ATG is the universal start codon in E. coli expression systems like pET-28a;
    TAA/TAG/TGA are the three stop codons in the standard genetic code.
    Rewarding their presence at the expected boundaries reduces false positives
    when multiple ORFs achieve the same alignment score.
    """
    if strand != "+":
        return 0

    dna = plasmid_dna
    L = len(dna)
    score = 0

    def _slice_circular(s: int, e: int) -> str:
        if e <= L and s <= e:
            return dna[s:e]
        # wrap
        return dna[s:] + dna[: (e % L)]

    start_codon = _slice_circular(start_nt, start_nt + 3)
    if start_codon == "ATG":
        score += 5

    stop_codon = _slice_circular(end_nt_exclusive, end_nt_exclusive + 3)
    if stop_codon in {"TAA", "TAG", "TGA"}:
        score += 5

    return score


def _match_with_x_wildcard(hay: str, needle: str) -> List[int]:
    """
    Return all start indices i where hay[i:i+len(needle)] matches needle,
    allowing 'X' in either string to match any AA.

    Note: this is only used when exact .find() fails (rare for clean DNA).
    """
    n = len(needle)
    if n == 0 or len(hay) < n:
        return []
    out: List[int] = []
    for i in range(0, len(hay) - n + 1):
        ok = True
        for j in range(n):
            a = hay[i + j]
            b = needle[j]
            if a != b and a != "X" and b != "X":
                ok = False
                break
        if ok:
            out.append(i)
    return out


def _best_fuzzy_identity(hay: str, needle: str) -> Tuple[float, int]:
    """
    Compute best same-length window identity between hay and needle.
    Returns (best_identity, best_start_index).

    'X' in either string is treated as a wildcard match (common with ambiguous DNA).
    """
    n = len(needle)
    if n == 0 or len(hay) < n:
        return 0.0, -1

    best_id = 0.0
    best_i = -1
    for i in range(0, len(hay) - n + 1):
        matches = 0
        for j in range(n):
            a = hay[i + j]
            b = needle[j]
            if a == b or a == "X" or b == "X":
                matches += 1
        ident = matches / n
        if ident > best_id:
            best_id = ident
            best_i = i
            if best_id >= 0.9999:  # early exit (near-perfect)
                break
    return best_id, best_i


# ----------------------------
# Smith–Waterman local alignment
# ----------------------------

def _smith_waterman_local_fallback(q: str, t: str, *, match: int = 2, mismatch: int = -1, gap: int = -2) -> Dict[str, Any]:
    """
    Classic Smith–Waterman local alignment (O(n*m)).
    Returns dict with:
      score, q_start, q_end, t_start, t_end, aligned_q, aligned_t, matches, aligned_query_len
    """
    n = len(q)
    m = len(t)
    # DP score matrix (only keep full for traceback clarity; guarded by size limits upstream)
    H = [[0] * (m + 1) for _ in range(n + 1)]
    P = [[0] * (m + 1) for _ in range(n + 1)]  # 0 stop, 1 diag, 2 up, 3 left

    best_score = 0
    best_pos = (0, 0)

    for i in range(1, n + 1):
        qi = q[i - 1]
        for j in range(1, m + 1):
            tj = t[j - 1]

            # Treat X as wildcard match (keeps alignment robust with ambiguous translation)
            is_match = (qi == tj) or (qi == "X") or (tj == "X")
            s_diag = H[i - 1][j - 1] + (match if is_match else mismatch)
            s_up = H[i - 1][j] + gap
            s_left = H[i][j - 1] + gap

            s = 0
            p = 0
            if s_diag >= s and s_diag > 0:
                s = s_diag
                p = 1
            if s_up >= s and s_up > 0:
                s = s_up
                p = 2
            if s_left >= s and s_left > 0:
                s = s_left
                p = 3

            H[i][j] = s
            P[i][j] = p

            if s > best_score:
                best_score = s
                best_pos = (i, j)

    # Traceback from best_pos
    i, j = best_pos
    aligned_q: List[str] = []
    aligned_t: List[str] = []
    matches = 0
    aligned_query_len = 0  # number of query residues aligned (non-gap in query)

    q_end = i
    t_end = j

    while i > 0 and j > 0 and P[i][j] != 0:
        p = P[i][j]
        if p == 1:  # diag
            qi = q[i - 1]
            tj = t[j - 1]
            aligned_q.append(qi)
            aligned_t.append(tj)
            if qi != "-" and tj != "-":
                aligned_query_len += 1
                if qi == tj or qi == "X" or tj == "X":
                    matches += 1
            i -= 1
            j -= 1
        elif p == 2:  # up
            qi = q[i - 1]
            aligned_q.append(qi)
            aligned_t.append("-")
            aligned_query_len += 1
            i -= 1
        else:  # left
            aligned_q.append("-")
            aligned_t.append(t[j - 1])
            j -= 1

    q_start = i
    t_start = j

    aligned_q.reverse()
    aligned_t.reverse()

    return {
        "score": best_score,
        "q_start": q_start,
        "q_end": q_end,
        "t_start": t_start,
        "t_end": t_end,
        "aligned_q": "".join(aligned_q),
        "aligned_t": "".join(aligned_t),
        "matches": matches,
        "aligned_query_len": aligned_query_len,
    }


def _smith_waterman_local(q: str, t: str, *, match: int = 2, mismatch: int = -1, gap: int = -2) -> Dict[str, Any]:
    """
    Smith–Waterman local alignment using Biopython when available (PairwiseAligner).
    Falls back to the in-house implementation if Biopython is missing.
    """
    try:
        from Bio.Align import PairwiseAligner
    except Exception:
        return _smith_waterman_local_fallback(q, t, match=match, mismatch=mismatch, gap=gap)

    # Configure local alignment
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap
    aligner.extend_gap_score = gap

    alignments = aligner.align(q, t)
    if not alignments:
        return _smith_waterman_local_fallback(q, t, match=match, mismatch=mismatch, gap=gap)

    aln = alignments[0]
    # aln.aligned is a tuple of (query_blocks, target_blocks)
    q_blocks, t_blocks = aln.aligned
    if len(q_blocks) == 0 or len(t_blocks) == 0:
        return _smith_waterman_local_fallback(q, t, match=match, mismatch=mismatch, gap=gap)

    q_start = int(q_blocks[0][0])
    q_end = int(q_blocks[-1][1])
    t_start = int(t_blocks[0][0])
    t_end = int(t_blocks[-1][1])

    aligned_query_len = 0
    matches = 0
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        aligned_query_len += int(qe - qs)
        for i in range(qe - qs):
            aq = q[qs + i]
            at = t[ts + i]
            if aq == at or aq == "X" or at == "X":
                matches += 1

    return {
        "score": int(aln.score),
        "q_start": q_start,
        "q_end": q_end,
        "t_start": t_start,
        "t_end": t_end,
        "aligned_q": "",
        "aligned_t": "",
        "matches": matches,
        "aligned_query_len": aligned_query_len,
    }


def _plausibility_warnings(plasmid_dna: str, strand: Strand, start_nt: int, end_nt_exclusive: int) -> List[str]:
    """
    Warnings-only checks to help catch "wrong region" calls.
    We intentionally do NOT fail valid engineered constructs; we only report.
    """
    warnings: List[str] = []
    if strand != "+":
        return warnings  # keep simple; '-' strand plausibility is less interpretable here

    dna = plasmid_dna
    L = len(dna)

    def _slice_circular(s: int, e: int) -> str:
        if e <= L and s <= e:
            return dna[s:e]
        return dna[s:] + dna[: (e % L)]

    try:
        start_codon = _slice_circular(start_nt, start_nt + 3)
        if start_codon != "ATG":
            warnings.append(f"Start codon at called start is not ATG (found {start_codon}).")

        stop_codon = _slice_circular(end_nt_exclusive, end_nt_exclusive + 3)
        if stop_codon and stop_codon not in {"TAA", "TAG", "TGA"}:
            warnings.append(f"No canonical stop codon immediately after ORF (found {stop_codon}).")
    except Exception:
        warnings.append("Plausibility checks could not be computed (non-fatal).")

    return warnings


# ----------------------------
# Main public API
# ----------------------------

def find_wt_in_plasmid(
    plasmid_dna: str,
    wt_protein: str,
    *,
    min_wt_len: int = 30,
    fuzzy_fallback: bool = True,
    min_identity: float = 0.95,
    align_min_identity: float = 0.90,
    align_min_coverage: float = 0.95,
    codon_table: Optional[Dict[str, str]] = None,
    max_align_wt_len: int = 2000,
    max_align_plasmid_len: int = 200000,
    allow_slow_alignment: bool = False,
    enable_plausibility_warnings: bool = True,
) -> GeneCall:
    """
    Validate that a circular plasmid encodes the given WT protein.

    Returns GeneCall with match_type:
      - "exact": exact AA substring found in a frame translation (X wildcard allowed)
      - "fuzzy": best same-length window identity passes threshold (handles substitutions)
      - "alignment": Smith–Waterman passes identity+coverage (handles indels)
      - "none": no candidate met thresholds
    """
    # if no codon table is provided, default to the standard genetic code —
    # keeps behaviour consistent between tests, CLI usage, and web routes
    if codon_table is None:
        codon_table = CODON_TABLE

    dna = _clean_dna(plasmid_dna)
    wt = _clean_protein(wt_protein)

    base_diag: Dict[str, Any] = {
        "params": {
            "min_wt_len": min_wt_len,
            "min_identity": min_identity,
            "align_min_identity": align_min_identity,
            "align_min_coverage": align_min_coverage,
            "max_align_wt_len": max_align_wt_len,
            "max_align_plasmid_len": max_align_plasmid_len,
            "allow_slow_alignment": allow_slow_alignment,
            "fuzzy_fallback": fuzzy_fallback,
        },
        "top_fuzzy_candidates": [],
        "top_alignment_candidates": [],
    }

    if len(wt) < min_wt_len:
        return GeneCall(
            is_valid=False,
            strand=None,
            frame=None,
            start_nt=None,
            end_nt_exclusive=None,
            wraps_origin=False,
            match_type="none",
            notes=f"WT protein too short for reliable validation (len={len(wt)} < {min_wt_len}).",
            diagnostics=base_diag,
        )

    if not dna:
        return GeneCall(
            is_valid=False,
            strand=None,
            frame=None,
            start_nt=None,
            end_nt_exclusive=None,
            wraps_origin=False,
            match_type="none",
            notes="Empty plasmid DNA provided.",
            diagnostics=base_diag,
        )

    L = len(dna)
    # doubling the plasmid sequence before translation exposes ORFs that cross
    # the origin — a gene cloned near the assembly start point would otherwise
    # be split across the boundary and missed by a linear search
    dna2 = dna + dna
    frames = translate_six_frames(dna2, codon_table=codon_table)

    # expected coding sequence length in nucleotides — used to map AA
    # match positions back to nucleotide coordinates on the original plasmid
    length_nt = len(wt) * 3

    # ----------------------------
    # A) Exact match (with X wildcard)
    # ----------------------------
    exact_candidates: List[_Candidate] = []

    for k, aa_seq in frames.items():
        strand, frame = _frame_key_to_strand_frame(k)

        # First try a fast exact .find() for clean sequences
        idxs: List[int] = []
        fast_i = aa_seq.find(wt)
        if fast_i != -1:
            idxs = [fast_i]
            x_note = ""
        else:
            # Fall back to wildcard scan ONLY if needed.
            idxs = _match_with_x_wildcard(aa_seq, wt)
            x_note = " X-wildcard enabled due to ambiguous translation." if idxs else ""

        for aa_i in idxs:
            if strand == "+":
                nt_start2 = frame + aa_i * 3
                nt_end2 = nt_start2 + length_nt
                start_nt, end_nt_exclusive, wraps_origin = _map_plus_nt_coords(L, nt_start2, length_nt)
            else:
                # Match is on rc(dna2); convert to plus-strand dna2 coordinates.
                nt_start_rev2 = frame + aa_i * 3
                start_nt, end_nt_exclusive, wraps_origin, nt_start2 = _map_minus_nt_coords(L, len(dna2), nt_start_rev2, length_nt)
                nt_end2 = nt_start2 + length_nt

            score = _orf_context_score(dna, strand, start_nt, end_nt_exclusive)

            exact_candidates.append(
                _Candidate(
                    strand=strand,
                    frame=frame,
                    frame_key=k,
                    aa_index=aa_i,
                    nt_start2=nt_start2,
                    nt_end2=nt_end2,
                    wraps_origin=wraps_origin,
                    score=score,
                    match_type="exact",
                    identity=1.0,
                    coverage=1.0,
                    notes=f"Exact WT match found in frame {strand}{frame}. start2={nt_start2}, end2={nt_end2} on doubled plasmid.{x_note}",
                )
            )

    best: Optional[_Candidate] = None
    if exact_candidates:
        # Choose best by ORF-context score, then earliest start on original plasmid for determinism
        exact_candidates.sort(key=lambda c: (-c.score, c.nt_start2))
        best = exact_candidates[0]
        if len(exact_candidates) > 1:
            best.notes += f" Multiple hits found ({len(exact_candidates)}); selected best by ORF-context score={best.score}."

    # ----------------------------
    # B) Fuzzy window identity (substitutions only)
    # ----------------------------
    # fuzzy fallback handles variants with point mutations relative to WT —
    # a 95% identity threshold allows up to ~5% substitution while still
    # confirming the correct gene is present rather than a random sequence
    if best is None and fuzzy_fallback:
        fuzzy_best: Optional[_Candidate] = None
        for k, aa_seq in frames.items():
            strand, frame = _frame_key_to_strand_frame(k)
            best_id, best_i = _best_fuzzy_identity(aa_seq, wt)
            base_diag["top_fuzzy_candidates"].append(
                {
                    "frame": k,
                    "best_identity": round(best_id, 6),
                    "passed": bool(best_id >= min_identity),
                    "used_x_wildcard": ("X" in aa_seq) or ("X" in wt),
                }
            )
            if best_i == -1:
                continue
            if best_id >= min_identity:
                if strand == "+":
                    nt_start2 = frame + best_i * 3
                    nt_end2 = nt_start2 + length_nt
                    start_nt, end_nt_exclusive, wraps_origin = _map_plus_nt_coords(L, nt_start2, length_nt)
                else:
                    nt_start_rev2 = frame + best_i * 3
                    start_nt, end_nt_exclusive, wraps_origin, nt_start2 = _map_minus_nt_coords(L, len(dna2), nt_start_rev2, length_nt)
                    nt_end2 = nt_start2 + length_nt

                score = _orf_context_score(dna, strand, start_nt, end_nt_exclusive)

                cand = _Candidate(
                    strand=strand,
                    frame=frame,
                    frame_key=k,
                    aa_index=best_i,
                    nt_start2=nt_start2,
                    nt_end2=nt_end2,
                    wraps_origin=wraps_origin,
                    score=score,
                    match_type="fuzzy",
                    identity=best_id,
                    coverage=1.0,
                    notes=f"Fuzzy match accepted in frame {strand}{frame}: identity={best_id:.3f} >= {min_identity:.3f}. start2={nt_start2}, end2={nt_end2}.",
                )
                if (fuzzy_best is None) or (cand.identity > fuzzy_best.identity) or (cand.identity == fuzzy_best.identity and cand.score > fuzzy_best.score):
                    fuzzy_best = cand

        # keep only a few diagnostic summaries
        base_diag["top_fuzzy_candidates"] = sorted(
            base_diag["top_fuzzy_candidates"], key=lambda d: d["best_identity"], reverse=True
        )[:3]

        best = fuzzy_best

    # ----------------------------
    # C) Smith–Waterman (indels) — guarded by size
    # ----------------------------
    # Smith-Waterman is the last resort — it handles indels that fuzzy matching
    # cannot, but is O(n*m) in time and memory so it is guarded by size limits
    # to prevent the server from hanging on large inputs
    if best is None:
        do_align = allow_slow_alignment or (len(wt) <= max_align_wt_len and (len(dna2) <= max_align_plasmid_len))
        if do_align:
            align_best: Optional[_Candidate] = None
            for k, aa_seq in frames.items():
                strand, frame = _frame_key_to_strand_frame(k)

                sw = _smith_waterman_local(wt, aa_seq)
                aligned_q_len = sw["aligned_query_len"]
                if aligned_q_len <= 0:
                    continue

                # identity over aligned query residues rather than full WT length —
                # avoids penalising a correct alignment that doesn't cover terminal regions
                identity = sw["matches"] / aligned_q_len
                # coverage confirms that most of the WT sequence was aligned,
                # preventing a short high-identity match from being accepted as valid
                coverage = aligned_q_len / len(wt)

                base_diag["top_alignment_candidates"].append(
                    {
                        "frame": k,
                        "identity": round(identity, 6),
                        "coverage": round(coverage, 6),
                        "sw_score": sw["score"],
                        "passed": bool(identity >= align_min_identity and coverage >= align_min_coverage),
                    }
                )

                if identity < align_min_identity or coverage < align_min_coverage:
                    continue

                # Map t_start (AA index in translated frame) to nt coords on dna2.
                # For alignment we map using the *aligned query length* (approx coding span on plasmid).
                # This is acceptable for coursework; downstream mutation calling will normally re-align.
                aa_i = sw["t_start"]
                length_nt_aln = aligned_q_len * 3

                if strand == "+":
                    nt_start2 = frame + aa_i * 3
                    nt_end2 = nt_start2 + length_nt_aln
                    start_nt, end_nt_exclusive, wraps_origin = _map_plus_nt_coords(L, nt_start2, length_nt_aln)
                else:
                    nt_start_rev2 = frame + aa_i * 3
                    start_nt, end_nt_exclusive, wraps_origin, nt_start2 = _map_minus_nt_coords(L, len(dna2), nt_start_rev2, length_nt_aln)
                    nt_end2 = nt_start2 + length_nt_aln

                score = _orf_context_score(dna, strand, start_nt, end_nt_exclusive)

                cand = _Candidate(
                    strand=strand,
                    frame=frame,
                    frame_key=k,
                    aa_index=aa_i,
                    nt_start2=nt_start2,
                    nt_end2=nt_end2,
                    wraps_origin=wraps_origin,
                    score=score,
                    match_type="align",
                    identity=identity,
                    coverage=coverage,
                    notes=(
                        f"Alignment match accepted in frame {strand}{frame}: "
                        f"identity={identity:.3f} (>= {align_min_identity:.3f}), "
                        f"coverage={coverage:.3f} (>= {align_min_coverage:.3f})."
                    ),
                )
                if (align_best is None) or (cand.coverage > align_best.coverage) or (
                    cand.coverage == align_best.coverage and cand.identity > align_best.identity
                ):
                    align_best = cand

            base_diag["top_alignment_candidates"] = sorted(
                base_diag["top_alignment_candidates"], key=lambda d: (d["coverage"], d["identity"]), reverse=True
            )[:3]

            best = align_best
        else:
            # record why alignment was skipped
            base_diag["alignment_skipped"] = {
                "reason": "Size guard triggered. Set allow_slow_alignment=True to force alignment.",
                "wt_len": len(wt),
                "dna2_len": len(dna2),
                "max_align_wt_len": max_align_wt_len,
                "max_align_plasmid_len": max_align_plasmid_len,
            }

    # ----------------------------
    # Final: return result
    # ----------------------------
    if best is None:
        return GeneCall(
            is_valid=False,
            strand=None,
            frame=None,
            start_nt=None,
            end_nt_exclusive=None,
            wraps_origin=False,
            match_type="none",
            identity=None,
            coverage=None,
            notes=(
                f"No exact match found. "
                f"Fuzzy fallback did not reach threshold (min_identity={min_identity:.2f}). "
                f"No alignment match met thresholds (identity>={align_min_identity:.2f}, coverage>={align_min_coverage:.2f})."
            ),
            diagnostics=base_diag,
        )

    # Attach plausibility warnings as diagnostics on success (warnings only)
    diagnostics: Optional[Dict[str, Any]] = None
    if enable_plausibility_warnings:
        warns = _plausibility_warnings(dna, best.strand, _map_plus_nt_coords(L, best.nt_start2, 0)[0] if best.strand == "+" else (best.nt_start2 % L),  # start_nt already in candidate mapping
                                       # best end_nt_exclusive is computed below; re-derive from mapping to avoid drift
                                       _map_plus_nt_coords(L, best.nt_start2, (len(wt) * 3 if best.match_type != "alignment" else int(round(best.coverage * len(wt))) * 3))[1])
        # The mapping above is conservative; we’ll compute canonical output coords below and re-check warnings.
        # We keep warnings if any.
        if warns:
            diagnostics = {"plausibility_warnings": warns}

    # Compute canonical output coords from candidate:
    if best.match_type in {"exact", "fuzzy"}:
        length_nt_out = len(wt) * 3
    else:
        length_nt_out = int(round(best.coverage * len(wt))) * 3

    # start/end on original plasmid:
    start_nt, end_nt_exclusive, wraps_origin = _map_plus_nt_coords(L, best.nt_start2, length_nt_out)

    # Recompute warnings with final coords:
    if enable_plausibility_warnings:
        warns = _plausibility_warnings(dna, best.strand, start_nt, end_nt_exclusive)
        if warns:
            diagnostics = diagnostics or {}
            diagnostics["plausibility_warnings"] = warns

    return GeneCall(
        is_valid=True,
        strand=best.strand,
        frame=best.frame,
        start_nt=start_nt,
        end_nt_exclusive=end_nt_exclusive,
        wraps_origin=wraps_origin,
        match_type=best.match_type,
        identity=best.identity,
        coverage=best.coverage,
        notes=best.notes,
        diagnostics=diagnostics,
    )
