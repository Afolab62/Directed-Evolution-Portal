"use client";

import { useEffect, useState } from "react";
import dynamic from "next/dynamic";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Mountain, Loader2, AlertCircle, RefreshCw } from "lucide-react";
import type { VariantData } from "@/lib/types";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false });

interface ActivityLandscapeProps {
  variants: VariantData[];
  experimentId: string;
}

type Method = "pca" | "tsne" | "umap";

interface LandscapeResponse {
  success: boolean;
  error?: string;
  figure: {
    data: Plotly.Data[];
    layout: Partial<Plotly.Layout>;
    frames?: Plotly.Frame[];
  };
  variant_count: number;
  method: Method;
}

export function ActivityLandscape({
  variants,
  experimentId,
}: ActivityLandscapeProps) {
  const [response, setResponse] = useState<LandscapeResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [method, setMethod] = useState<Method>("pca");

  // Track Plotly's internal figure state so slider drags, play/stop presses,
  // and legend clicks are not clobbered by React re-renders.
  const [figData, setFigData] = useState<Plotly.Data[]>([]);
  const [figLayout, setFigLayout] = useState<Partial<Plotly.Layout>>({});
  const [figFrames, setFigFrames] = useState<Plotly.Frame[]>([]);

  // Whenever a fresh landscape is fetched, seed the tracked state.
  useEffect(() => {
    if (response?.figure) {
      setFigData(response.figure.data);
      setFigLayout({ ...response.figure.layout, autosize: true });
      setFigFrames(response.figure.frames ?? []);
    }
  }, [response]);

  const fetchLandscape = async (selectedMethod: Method) => {
    setLoading(true);
    setError(null);
    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${experimentId}/landscape?method=${selectedMethod}`,
        { credentials: "include" },
      );
      const json: LandscapeResponse = await res.json();
      if (!json.success)
        throw new Error(json.error ?? "Failed to load landscape");
      setResponse(json);
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : "Unknown error");
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    if (experimentId && variants.length >= 3) fetchLandscape(method);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [experimentId]);

  const handleMethodChange = (m: Method) => {
    setMethod(m);
    fetchLandscape(m);
  };

  if (variants.length === 0) {
    return (
      <Card>
        <CardContent className="py-12 text-center">
          <Mountain className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
          <p className="text-muted-foreground">
            No data available for landscape visualization
          </p>
        </CardContent>
      </Card>
    );
  }

  const analysedCount = variants.filter(
    (v) =>
      v.proteinSequence ||
      (v.mutationCount !== undefined && v.mutationCount > 0),
  ).length;

  if (analysedCount < 3) {
    return (
      <Card>
        <CardContent className="py-12 text-center">
          <AlertCircle className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
          <p className="text-muted-foreground font-medium">
            Sequence analysis required
          </p>
          <p className="text-sm text-muted-foreground mt-1">
            {analysedCount} of {variants.length} variants have been analysed.
            Run sequence analysis first.
          </p>
        </CardContent>
      </Card>
    );
  }

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <div>
            <CardTitle className="flex items-center gap-2">
              <Mountain className="h-5 w-5" />
              Activity Landscape
            </CardTitle>
            <CardDescription>
              3D topographic surface — peaks represent highly active variants.
              Use the play buttons and slider to animate across generations.
            </CardDescription>
          </div>
          <div className="flex items-center gap-2">
            {(["pca", "tsne", "umap"] as Method[]).map((m) => (
              <Button
                key={m}
                size="sm"
                variant={method === m ? "default" : "outline"}
                onClick={() => handleMethodChange(m)}
                disabled={loading}
              >
                {m.toUpperCase()}
              </Button>
            ))}
            {response && (
              <Button
                size="sm"
                variant="ghost"
                onClick={() => fetchLandscape(method)}
                disabled={loading}
              >
                <RefreshCw className="h-4 w-4" />
              </Button>
            )}
          </div>
        </div>
      </CardHeader>

      <CardContent>
        {loading && (
          <div className="flex items-center justify-center h-64 text-muted-foreground">
            <Loader2 className="h-6 w-6 animate-spin mr-2" />
            Computing {method.toUpperCase()} landscape…
          </div>
        )}

        {error && !loading && (
          <div className="flex flex-col items-center justify-center h-40 gap-3">
            <AlertCircle className="h-8 w-8 text-destructive" />
            <p className="text-sm text-destructive text-center max-w-md">
              {error}
            </p>
            <Button
              size="sm"
              variant="outline"
              onClick={() => fetchLandscape(method)}
            >
              Retry
            </Button>
          </div>
        )}

        {response && !loading && (
          <>
            <p className="text-xs text-muted-foreground mb-2">
              {response.variant_count} variants ·{" "}
              {response.method.toUpperCase()} projection · use ▶ buttons or
              slider below to animate generations
            </p>
            <Plot
              data={figData}
              layout={figLayout}
              frames={figFrames}
              onInitialized={(fig) => {
                setFigData(fig.data as Plotly.Data[]);
                setFigLayout(fig.layout as Partial<Plotly.Layout>);
              }}
              onUpdate={(fig) => {
                setFigData(fig.data as Plotly.Data[]);
                setFigLayout(fig.layout as Partial<Plotly.Layout>);
              }}
              config={{ responsive: true, displayModeBar: true }}
              style={{ width: "100%" }}
              useResizeHandler
            />
          </>
        )}

        {!response && !loading && !error && (
          <div className="flex items-center justify-center h-40">
            <p className="text-muted-foreground text-sm">
              Click a method above to compute the landscape
            </p>
          </div>
        )}
      </CardContent>
    </Card>
  );
}
