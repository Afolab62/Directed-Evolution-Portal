"use client";

import { useEffect, useMemo, useRef, useState } from "react";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Skeleton } from "@/components/ui/skeleton";
import { BarChart3 } from "lucide-react";
import type { VariantData } from "@/lib/types";

const BACKEND = process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000";

interface ActivityDistributionChartProps {
  variants: VariantData[];
  experimentId: string;
}

export function ActivityDistributionChart({
  variants,
  experimentId,
}: ActivityDistributionChartProps) {
  const [imgSrc, setImgSrc] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const prevUrl = useRef<string | null>(null);

  useEffect(() => {
    if (!experimentId) return;
    let active = true;
    setLoading(true);
    setError(null);

    fetch(
      `${BACKEND}/api/experiments/${experimentId}/plots/activity-distribution`,
      { credentials: "include" },
    )
      .then((res) => {
        if (!res.ok)
          return res.json().then((j) => {
            throw new Error(j.error || `Server error ${res.status}`);
          });
        return res.blob();
      })
      .then((blob) => {
        if (!active) return;
        // Revoke previous blob URL to avoid memory leak
        if (prevUrl.current) URL.revokeObjectURL(prevUrl.current);
        const url = URL.createObjectURL(blob);
        prevUrl.current = url;
        setImgSrc(url);
      })
      .catch((err) => {
        if (active) setError(err.message);
      })
      .finally(() => {
        if (active) setLoading(false);
      });

    return () => {
      active = false;
    };
  }, [experimentId]);

  // Summary stats computed from variants prop (unchanged)
  const generationStats = useMemo(() => {
    const generations = [...new Set(variants.map((v) => v.generation))].sort(
      (a, b) => a - b,
    );

    return generations.map((gen) => {
      const scores = variants
        .filter((v) => v.generation === gen)
        .map((v) => v.activityScore)
        .filter((s) => s != null && isFinite(s))
        .sort((a, b) => a - b);

      if (scores.length === 0)
        return {
          generation: gen,
          count: 0,
          meanActivity: 0,
          medianActivity: 0,
          minActivity: 0,
          maxActivity: 0,
          stdDev: 0,
        };

      const mean = scores.reduce((s, v) => s + v, 0) / scores.length;
      const mid = Math.floor(scores.length / 2);
      const median =
        scores.length % 2 === 0
          ? (scores[mid - 1] + scores[mid]) / 2
          : scores[mid];
      const stdDev = Math.sqrt(
        scores.reduce((s, v) => s + Math.pow(v - mean, 2), 0) / scores.length,
      );

      return {
        generation: gen,
        count: scores.length,
        meanActivity: mean,
        medianActivity: median,
        minActivity: scores[0],
        maxActivity: scores[scores.length - 1],
        stdDev,
      };
    });
  }, [variants]);

  if (generationStats.length === 0) {
    return (
      <Card>
        <CardContent className="py-12 text-center">
          <BarChart3 className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
          <p className="text-muted-foreground">No generation data available</p>
        </CardContent>
      </Card>
    );
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle className="flex items-center gap-2">
          <BarChart3 className="h-5 w-5" />
          Activity Score Distribution by Generation
        </CardTitle>
        <CardDescription>
          Violin plot showing the full score distribution, mean, and median per
          generation
        </CardDescription>
      </CardHeader>
      <CardContent>
        {loading && <Skeleton className="w-full h-95 rounded-md" />}
        {error && !loading && (
          <div className="w-full h-95 flex items-center justify-center text-muted-foreground text-sm">
            Failed to load plot: {error}
          </div>
        )}
        {imgSrc && !loading && (
          <img
            src={imgSrc}
            alt="Activity score distribution by generation"
            className="w-full max-h-[500px] object-contain rounded-lg border border-border"
          />
        )}

        {/* Generation summary table */}
        <div className="mt-4 overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b">
                <th className="text-left py-2 px-2 font-medium">Generation</th>
                <th className="text-left py-2 px-2 font-medium">Variants</th>
                <th className="text-left py-2 px-2 font-medium">Mean</th>
                <th className="text-left py-2 px-2 font-medium">Median</th>
                <th className="text-left py-2 px-2 font-medium">Min</th>
                <th className="text-left py-2 px-2 font-medium">Max</th>
                <th className="text-left py-2 px-2 font-medium">Std Dev</th>
              </tr>
            </thead>
            <tbody>
              {generationStats.map((stat) => (
                <tr key={stat.generation} className="border-b">
                  <td className="py-2 px-2 font-medium">
                    Gen {stat.generation}
                  </td>
                  <td className="py-2 px-2">{stat.count}</td>
                  <td className="py-2 px-2 font-mono">
                    {stat.meanActivity.toFixed(3)}
                  </td>
                  <td className="py-2 px-2 font-mono">
                    {stat.medianActivity.toFixed(3)}
                  </td>
                  <td className="py-2 px-2 font-mono text-muted-foreground">
                    {stat.minActivity.toFixed(3)}
                  </td>
                  <td className="py-2 px-2 font-mono text-green-400">
                    {stat.maxActivity.toFixed(3)}
                  </td>
                  <td className="py-2 px-2 font-mono text-muted-foreground">
                    {stat.stdDev.toFixed(3)}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </CardContent>
    </Card>
  );
}
