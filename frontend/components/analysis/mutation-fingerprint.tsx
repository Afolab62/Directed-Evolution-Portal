"use client";

import { useState, useEffect } from "react";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Tabs, TabsList, TabsTrigger, TabsContent } from "@/components/ui/tabs";
import {
  Loader2,
  Dna,
  AlertCircle,
  Download,
  BarChart2,
  Box,
  Crosshair,
  X,
} from "lucide-react";
import type { VariantData } from "@/lib/types";

const BACKEND = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";

// ── Props ─────────────────────────────────────────────────────────────────

interface MutationFingerprintProps {
  variants: VariantData[];
  experimentId: string;
  selectedVariantIndex: number | null;
  onSelectVariant: (index: number) => void;
  wtSequence?: string;
}

// ── Iframe figure component ───────────────────────────────────────────────
// Points an <iframe src> directly at the backend ?format=html endpoint.
// The browser fetches and renders it natively — no JS string wrangling,
// no srcDoc size limits, CDN scripts load freely in the iframe's own context.

function PlotlyIframe({
  url,
  height,
  label,
}: {
  url: string;
  height: number;
  label: string;
}) {
  const [loaded, setLoaded] = useState(false);
  const [error, setError] = useState(false);

  // Reset state whenever the URL changes (new variant selected)
  useEffect(() => {
    setLoaded(false);
    setError(false);
  }, [url]);

  return (
    <div style={{ position: "relative", width: "100%", height }}>
      {/* Spinner shown until iframe fires onLoad */}
      {!loaded && !error && (
        <div
          className="absolute inset-0 flex items-center justify-center gap-3"
          style={{ zIndex: 1 }}
        >
          <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
          <span className="text-muted-foreground text-sm">
            Loading {label}…
          </span>
        </div>
      )}
      {error && (
        <div
          className="absolute inset-0 flex flex-col items-center justify-center gap-2"
          style={{ zIndex: 1 }}
        >
          <AlertCircle className="h-8 w-8 text-destructive" />
          <p className="text-sm text-destructive">Failed to load {label}</p>
        </div>
      )}
      <iframe
        key={url}
        src={url}
        width="100%"
        height={height}
        style={{ border: "none", display: "block", opacity: loaded ? 1 : 0 }}
        title={label}
        onLoad={() => setLoaded(true)}
        onError={() => setError(true)}
      />
    </div>
  );
}

// ── Main component ────────────────────────────────────────────────────────

export function MutationFingerprint({
  variants,
  experimentId,
  selectedVariantIndex,
  onSelectVariant,
  wtSequence,
}: MutationFingerprintProps) {
  const [activeTab, setActiveTab] = useState<"linear" | "3d">("linear");
  const [fp3dRequested, setFp3dRequested] = useState(false);
  const [highlightedPosition, setHighlightedPosition] = useState<number | null>(
    null,
  );

  // Reset when variant changes
  useEffect(() => {
    setFp3dRequested(false);
    setActiveTab("linear");
    setHighlightedPosition(null);
  }, [selectedVariantIndex]);

  // Listen for clicks forwarded from the linear fingerprint iframe.
  // When the user clicks a mutation triangle, the iframe fires:
  //   window.parent.postMessage({ type: 'dem_mutation_click', position: N }, '*')
  // We catch it here, record the position, and switch to the 3D tab with the
  // corresponding residue highlighted as a gold sphere.
  useEffect(() => {
    const handle = (e: MessageEvent) => {
      if (!e.data || e.data.type !== "dem_mutation_click") return;
      const pos = Number(e.data.position);
      if (!Number.isFinite(pos)) return;
      setHighlightedPosition(pos);
      setFp3dRequested(true);
      setActiveTab("3d");
    };
    window.addEventListener("message", handle);
    return () => window.removeEventListener("message", handle);
  }, []);

  const selectedVariant =
    variants.find((v) => v.plasmidVariantIndex === selectedVariantIndex) ??
    null;

  // Iframe URLs — backend returns self-contained Plotly HTML for ?format=html
  const linearHtmlUrl = selectedVariant
    ? `${BACKEND}/api/experiments/${experimentId}/fingerprint_linear/${selectedVariant.id}?format=html`
    : null;
  const fp3dHtmlUrl = selectedVariant
    ? `${BACKEND}/api/experiments/${experimentId}/fingerprint3d/${selectedVariant.id}?format=html${
        highlightedPosition != null ? `&highlight=${highlightedPosition}` : ""
      }`
    : null;

  const handleTabChange = (tab: string) => {
    setActiveTab(tab as "linear" | "3d");
    if (tab === "3d") setFp3dRequested(true);
  };

  // Summary counts derived from already-loaded VariantData — no extra fetch needed
  const mutationCount =
    selectedVariant?.mutationCount ??
    selectedVariant?.mutations?.length ??
    null;

  const synonymousCount =
    selectedVariant?.mutations?.filter((m) => m.type === "synonymous").length ??
    null;
  const nonSynonymousCount =
    selectedVariant?.mutations?.filter((m) => m.type === "non-synonymous")
      .length ?? null;

  // ── Render ───────────────────────────────────────────────────────────────
  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Dna className="h-5 w-5" />
          Mutation Fingerprint
        </CardTitle>
        <CardDescription>
          Linear view shows mutations along the sequence backbone by generation.
          3D view maps them onto the AlphaFold structure (pLDDT-coloured
          backbone).
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Variant selector */}
        <div>
          <label className="text-sm font-medium mb-2 block">
            Select Variant
          </label>
          <Select
            value={selectedVariantIndex?.toString() ?? ""}
            onValueChange={(v) => onSelectVariant(parseFloat(v))}
          >
            <SelectTrigger>
              <SelectValue placeholder="Choose a variant…" />
            </SelectTrigger>
            <SelectContent>
              {variants.map((v, idx) => (
                <SelectItem key={v.id} value={v.plasmidVariantIndex.toString()}>
                  #{idx + 1} — Variant {v.plasmidVariantIndex} | Gen{" "}
                  {v.generation} | Activity {(v.activityScore ?? 0).toFixed(3)}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>

        {!selectedVariant ? (
          <div className="flex flex-col items-center justify-center py-12 text-center">
            <Dna className="h-10 w-10 text-muted-foreground/40 mb-3" />
            <p className="text-muted-foreground">
              Select a variant above to view its mutation fingerprint
            </p>
          </div>
        ) : (
          <>
            {/* Summary badges — built from already-loaded VariantData, zero extra fetches */}
            <div className="flex flex-wrap items-center gap-3 text-sm">
              <span className="text-muted-foreground">
                Variant{" "}
                <span className="font-mono font-medium text-foreground">
                  {selectedVariant.plasmidVariantIndex}
                </span>
              </span>
              <Badge variant="outline">Gen {selectedVariant.generation}</Badge>
              <Badge variant="outline">
                Activity {(selectedVariant.activityScore ?? 0).toFixed(3)}
              </Badge>
              {mutationCount !== null && (
                <Badge variant="outline">
                  {mutationCount} mutation{mutationCount !== 1 ? "s" : ""}
                </Badge>
              )}
              {nonSynonymousCount !== null && (
                <Badge variant="outline">
                  {nonSynonymousCount} non-synonymous
                </Badge>
              )}
              {synonymousCount !== null && (
                <Badge variant="outline">{synonymousCount} synonymous</Badge>
              )}
              <Button
                size="sm"
                variant="outline"
                className="ml-auto gap-1.5"
                asChild
              >
                <a
                  href={`${BACKEND}/api/experiments/${experimentId}/mutations/export`}
                  download
                >
                  <Download className="h-3.5 w-3.5" />
                  Export CSV
                </a>
              </Button>
            </div>

            {/* Linear / 3D tabs */}
            <Tabs value={activeTab} onValueChange={handleTabChange}>
              <TabsList>
                <TabsTrigger value="linear" className="gap-1.5">
                  <BarChart2 className="h-4 w-4" />
                  Linear
                </TabsTrigger>
                <TabsTrigger value="3d" className="gap-1.5">
                  <Box className="h-4 w-4" />
                  3D Structure
                </TabsTrigger>
              </TabsList>

              {/* ── Linear tab ── */}
              <TabsContent value="linear" className="mt-4">
                {linearHtmlUrl && (
                  <PlotlyIframe
                    url={linearHtmlUrl}
                    height={560}
                    label="linear fingerprint"
                  />
                )}
              </TabsContent>

              {/* ── 3D Structure tab ── */}
              <TabsContent value="3d" className="mt-4 space-y-3">
                {!fp3dRequested ? (
                  <div className="flex flex-col items-center justify-center py-16 gap-3">
                    <Box className="h-10 w-10 text-muted-foreground/40" />
                    <p className="text-sm text-muted-foreground">
                      Click a mutation triangle in the Linear tab to jump
                      directly to it in 3D, or switch to this tab to load the
                      full structure.
                    </p>
                  </div>
                ) : fp3dHtmlUrl ? (
                  <>
                    {/* Highlight indicator — shown when user clicked a triangle in the linear view */}
                    {highlightedPosition != null ? (
                      <div className="flex items-center gap-2 px-3 py-2 rounded-md bg-yellow-50 border border-yellow-200 text-yellow-800 dark:bg-yellow-950 dark:border-yellow-800 dark:text-yellow-300">
                        <Crosshair className="h-4 w-4 shrink-0 text-yellow-600" />
                        <span className="text-sm flex-1">
                          Highlighting residue{" "}
                          <span className="font-mono font-semibold">
                            {highlightedPosition}
                          </span>{" "}
                          — gold sphere in the 3D view
                        </span>
                        <button
                          onClick={() => setHighlightedPosition(null)}
                          className="shrink-0 text-yellow-600 hover:text-yellow-800 dark:text-yellow-400"
                        >
                          <X className="h-4 w-4" />
                        </button>
                      </div>
                    ) : null}
                    <PlotlyIframe
                      url={fp3dHtmlUrl}
                      height={800}
                      label="3D structure fingerprint"
                    />
                  </>
                ) : null}
              </TabsContent>
            </Tabs>
          </>
        )}
      </CardContent>
    </Card>
  );
}
