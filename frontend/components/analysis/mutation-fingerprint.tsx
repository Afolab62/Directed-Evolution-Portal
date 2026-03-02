"use client";

import { useState, useEffect, useMemo } from "react";
import dynamic from "next/dynamic";
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
import { Loader2, Dna, AlertCircle } from "lucide-react";
import type { VariantData } from "@/lib/types";

// Plotly is browser-only â€” dynamic import
const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

// â”€â”€ Types â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
interface Fingerprint3dResponse {
  success: boolean;
  error?: string;
  // Plotly figure JSON returned by fig.to_json()
  figure: {
    data: Plotly.Data[];
    layout: Partial<Plotly.Layout>;
  };
  variantId: string;
  plasmidVariantIndex: number;
  generation: number;
  activityScore: number | null;
  numMutations: number;
  numSynonymous: number;
  numNonSynonymous: number;
  structureStatus: string;
  lineageLength: number;
}

// â”€â”€ Props â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
interface MutationFingerprintProps {
  variants: VariantData[];
  experimentId: string;
  selectedVariantIndex: number | null;
  onSelectVariant: (index: number) => void;
  wtSequence?: string; // kept for backward compat
}

export function MutationFingerprint({
  variants,
  experimentId,
  selectedVariantIndex,
  onSelectVariant,
}: MutationFingerprintProps) {
  const [fpData, setFpData] = useState<Fingerprint3dResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const selectedVariant = useMemo(
    () => variants.find((v) => v.plasmidVariantIndex === selectedVariantIndex) ?? null,
    [variants, selectedVariantIndex],
  );

  // Fetch the Plotly figure from the new 3D endpoint whenever selection changes
  useEffect(() => {
    if (!selectedVariant || !experimentId) {
      setFpData(null);
      return;
    }
    let cancelled = false;
    setLoading(true);
    setError(null);

    fetch(
      `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${experimentId}/fingerprint3d/${selectedVariant.id}`,
      { credentials: "include" },
    )
      .then((r) => r.json())
      .then((data: Fingerprint3dResponse) => {
        if (cancelled) return;
        if (data.success) setFpData(data);
        else setError(data.error ?? "Failed to load fingerprint");
      })
      .catch(() => {
        if (!cancelled) setError("Network error");
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });

    return () => {
      cancelled = true;
    };
  }, [selectedVariant?.id, experimentId]);

  // â”€â”€ Responsive layout override â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  // The backend produces a fixed height=860 layout; keep that but make it
  // fill the card width by overriding layout.autosize.
  const plotLayout = useMemo((): Partial<Plotly.Layout> => {
    if (!fpData?.figure?.layout) return {};
    return {
      ...fpData.figure.layout,
      autosize: true,
    };
  }, [fpData]);

  // â”€â”€ Render â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <Dna className="h-5 w-5" />
          Mutation Fingerprint
        </CardTitle>
        <CardDescription>
          Select a variant to view mutations mapped onto the 3D protein structure
          (AlphaFold backbone coloured by pLDDT confidence)
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Variant selector */}
        <div>
          <label className="text-sm font-medium mb-2 block">Select Variant</label>
          <Select
            value={selectedVariantIndex?.toString() ?? ""}
            onValueChange={(v) => onSelectVariant(parseFloat(v))}
          >
            <SelectTrigger>
              <SelectValue placeholder="Choose a top-performing variantâ€¦" />
            </SelectTrigger>
            <SelectContent>
              {variants.map((v, idx) => (
                <SelectItem key={v.id} value={v.plasmidVariantIndex.toString()}>
                  #{idx + 1} â€” Variant {v.plasmidVariantIndex} | Gen {v.generation} | Activity{" "}
                  {(v.activityScore ?? 0).toFixed(3)}
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
        ) : loading ? (
          <div className="flex flex-col items-center justify-center py-12 gap-3">
            <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
            <span className="text-muted-foreground text-sm">
              Fetching 3D structure…
            </span>
          </div>
        ) : error ? (
          <div className="flex flex-col items-center py-10 text-center">
            <AlertCircle className="h-8 w-8 text-destructive mb-2" />
            <p className="text-sm text-destructive">{error}</p>
          </div>
        ) : fpData && fpData.numMutations === 0 ? (
          <div className="text-center py-10">
            <p className="text-muted-foreground font-medium">
              Wild-type â€” no mutations in lineage
            </p>
            <p className="text-sm text-muted-foreground mt-1">
              This variant and its ancestors carry no changes relative to the WT
              reference.
            </p>
          </div>
        ) : fpData ? (
          <>
            {/* Summary badges */}
            <div className="flex flex-wrap gap-3 text-sm">
              <span className="text-muted-foreground">
                Variant{" "}
                <span className="font-mono font-medium text-foreground">
                  {fpData.plasmidVariantIndex}
                </span>
              </span>
              <Badge variant="outline">Gen {fpData.generation}</Badge>
              <Badge variant="outline">
                Activity {(fpData.activityScore ?? 0).toFixed(3)}
              </Badge>
              <Badge variant="outline">
                {fpData.numMutations} mutation
                {fpData.numMutations !== 1 ? "s" : ""}
              </Badge>
              <Badge variant="outline">
                {fpData.numNonSynonymous} non-synonymous
              </Badge>
              <Badge variant="outline">{fpData.numSynonymous} synonymous</Badge>
              <Badge variant="outline">
                {fpData.lineageLength} step
                {fpData.lineageLength !== 1 ? "s" : ""} in lineage
              </Badge>
              {fpData.structureStatus && (
                <Badge variant="secondary" className="text-xs">
                  {fpData.structureStatus}
                </Badge>
              )}
            </div>

            {/* Plotly circular fingerprint + optional 3D panel */}
            <div className="w-full overflow-x-auto">
              <Plot
                data={fpData.figure.data}
                layout={plotLayout}
                config={{ displayModeBar: true, responsive: true }}
                style={{ width: "100%" }}
                useResizeHandler
              />
            </div>
          </>
        ) : null}
      </CardContent>
    </Card>
  );
}
