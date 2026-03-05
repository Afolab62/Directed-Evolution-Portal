"use client";

import { useEffect, useState, Suspense } from "react";
import { useSearchParams } from "next/navigation";
import dynamic from "next/dynamic";
import { useAuth } from "@/lib/auth-context";
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
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Badge } from "@/components/ui/badge";
import { BarChart3, Loader2, TrendingUp, Dna, AlertCircle } from "lucide-react";
import { Skeleton } from "@/components/ui/skeleton";
import type { Experiment, VariantData } from "@/lib/types";
import Loading from "./loading";

// Dynamically import all heavy analysis components so Next.js code-splits
// them into separate chunks. This prevents react-plotly.js (~3 MB) and
// recharts from being parsed as part of the main page bundle, cutting the
// first-compile time from ~20 s down to ~3 s.
const TopPerformersTable = dynamic(
  () =>
    import("@/components/analysis/top-performers-table").then((m) => ({
      default: m.TopPerformersTable,
    })),
  { ssr: false, loading: () => <Skeleton className="h-48 w-full" /> },
);
const ActivityDistributionChart = dynamic(
  () =>
    import("@/components/analysis/activity-distribution-chart").then((m) => ({
      default: m.ActivityDistributionChart,
    })),
  { ssr: false, loading: () => <Skeleton className="h-64 w-full" /> },
);
const MutationFingerprint = dynamic(
  () =>
    import("@/components/analysis/mutation-fingerprint").then((m) => ({
      default: m.MutationFingerprint,
    })),
  { ssr: false, loading: () => <Skeleton className="h-96 w-full" /> },
);
const ActivityLandscape = dynamic(
  () =>
    import("@/components/analysis/activity-landscape").then((m) => ({
      default: m.ActivityLandscape,
    })),
  { ssr: false, loading: () => <Skeleton className="h-96 w-full" /> },
);

function AnalysisContent() {
  const searchParams = useSearchParams();
  const { user } = useAuth();

  const [experiments, setExperiments] = useState<Experiment[]>([]);
  const [selectedExpId, setSelectedExpId] = useState<string>(
    searchParams.get("experiment") || "",
  );
  const [variants, setVariants] = useState<VariantData[]>([]);
  const [topPerformers, setTopPerformers] = useState<VariantData[]>([]);
  const [selectedExperiment, setSelectedExperiment] =
    useState<Experiment | null>(null);
  const [isLoading, setIsLoading] = useState(true);
  const [isLoadingTopPerformers, setIsLoadingTopPerformers] = useState(false);
  const [selectedVariantIndex, setSelectedVariantIndex] = useState<
    number | null
  >(null);
  const [activeTab, setActiveTab] = useState("overview");
  const [loadError, setLoadError] = useState<string | null>(null);
  const [variantError, setVariantError] = useState<string | null>(null);

  // Load experiments list
  useEffect(() => {
    async function loadExperiments() {
      if (!user) return;
      setLoadError(null);

      try {
        const res = await fetch(
          `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments?userId=${user.id}`,
          {
            credentials: "include",
          },
        );
        if (!res.ok) throw new Error(`Server error ${res.status}`);
        const data = await res.json();
        if (data.success) {
          const validExps = data.experiments.filter(
            (e: Experiment) => e.validationStatus === "valid",
          );
          setExperiments(validExps);

          // Auto-select if only one or from URL param
          if (validExps.length === 1) {
            setSelectedExpId(validExps[0].id);
          } else if (searchParams.get("experiment")) {
            setSelectedExpId(searchParams.get("experiment") || "");
          }
        }
      } catch (err) {
        console.error("Failed to load experiments:", err);
        setLoadError(
          err instanceof Error
            ? err.message
            : "Failed to load experiments. Please refresh the page.",
        );
      } finally {
        setIsLoading(false);
      }
    }

    loadExperiments();
  }, [user, searchParams]);

  // Load variants when experiment selected
  useEffect(() => {
    async function loadVariants() {
      if (!selectedExpId) {
        setVariants([]);
        setSelectedExperiment(null);
        setVariantError(null);
        return;
      }
      setVariantError(null);
      setSelectedVariantIndex(null);
      setActiveTab("overview");

      try {
        const res = await fetch(
          `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${selectedExpId}`,
          {
            credentials: "include",
          },
        );
        if (!res.ok) throw new Error(`Server error ${res.status}`);
        const data = await res.json();

        if (data.success) {
          setSelectedExperiment(data.experiment);
          setVariants(data.variants || []);
        } else {
          throw new Error(data.error || "Failed to load experiment data");
        }
      } catch (err) {
        console.error("Failed to load variants:", err);
        setVariantError(
          err instanceof Error ? err.message : "Failed to load variant data.",
        );
      }
    }

    loadVariants();
  }, [selectedExpId]);

  // Load top performers with mutations
  useEffect(() => {
    async function loadTopPerformers() {
      if (!selectedExpId) {
        setTopPerformers([]);
        return;
      }

      setIsLoadingTopPerformers(true);
      try {
        const res = await fetch(
          `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${selectedExpId}/top-performers?limit=10&include_mutations=true`,
          {
            credentials: "include",
          },
        );
        const data = await res.json();

        if (data.success) {
          setTopPerformers(data.topPerformers || []);
        }
      } catch (err) {
        console.error("Failed to load top performers:", err);
      } finally {
        setIsLoadingTopPerformers(false);
      }
    }

    loadTopPerformers();
  }, [selectedExpId]);

  if (isLoading) {
    return null;
  }

  if (loadError) {
    return (
      <Card>
        <CardContent className="py-12 text-center">
          <AlertCircle className="h-12 w-12 mx-auto text-destructive/70 mb-3" />
          <p className="text-destructive font-medium">
            Failed to load experiments
          </p>
          <p className="text-muted-foreground text-sm mt-1">{loadError}</p>
        </CardContent>
      </Card>
    );
  }

  const passedVariants = variants.filter((v) => v.qcStatus === "passed");

  return (
    <div className="space-y-6">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <h1 className="text-2xl font-bold text-foreground">
            Analysis Dashboard
          </h1>
          <p className="text-muted-foreground mt-1">
            Visualise and analyse your directed evolution experiments
          </p>
        </div>

        <Select value={selectedExpId} onValueChange={setSelectedExpId}>
          <SelectTrigger className="w-full sm:w-[280px]">
            <SelectValue placeholder="Select experiment" />
          </SelectTrigger>
          <SelectContent>
            {experiments.map((exp) => (
              <SelectItem key={exp.id} value={exp.id}>
                {exp.name}
              </SelectItem>
            ))}
          </SelectContent>
        </Select>
      </div>

      {experiments.length === 0 ? (
        <Card>
          <CardContent className="py-12 text-center">
            <AlertCircle className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
            <p className="text-muted-foreground">
              No experiments with uploaded data found. Create an experiment and
              upload data to see analysis.
            </p>
          </CardContent>
        </Card>
      ) : !selectedExpId ? (
        <Card>
          <CardContent className="py-12 text-center">
            <BarChart3 className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
            <p className="text-muted-foreground">
              Select an experiment above to view analysis
            </p>
          </CardContent>
        </Card>
      ) : variantError ? (
        <Card>
          <CardContent className="py-12 text-center">
            <AlertCircle className="h-12 w-12 mx-auto text-destructive/70 mb-3" />
            <p className="text-destructive font-medium">
              Failed to load experiment data
            </p>
            <p className="text-muted-foreground text-sm mt-1">{variantError}</p>
          </CardContent>
        </Card>
      ) : variants.length === 0 ? (
        <Card>
          <CardContent className="py-12 text-center">
            <Dna className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
            <p className="text-muted-foreground">
              No variant data uploaded for this experiment. Upload data to see
              analysis.
            </p>
          </CardContent>
        </Card>
      ) : (
        <>
          // Summary stats //
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <Card>
              <CardContent className="pt-6">
                <p className="text-sm text-muted-foreground">Total Variants</p>
                <p className="text-2xl font-bold">{variants.length}</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-6">
                <p className="text-sm text-muted-foreground">Passed QC</p>
                <p className="text-2xl font-bold text-accent">
                  {passedVariants.length}
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-6">
                <p className="text-sm text-muted-foreground">Avg Activity</p>
                <p className="text-2xl font-bold">
                  {(
                    passedVariants.reduce((s, v) => s + v.activityScore, 0) /
                      passedVariants.length || 0
                  ).toFixed(2)}
                </p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-6">
                <div className="flex items-center gap-2">
                  <TrendingUp className="h-4 w-4 text-accent" />
                  <p className="text-sm text-muted-foreground">Best Activity</p>
                </div>
                <p className="text-2xl font-bold">
                  {Math.max(
                    ...passedVariants.map((v) => v.activityScore),
                    0,
                  ).toFixed(2)}
                </p>
              </CardContent>
            </Card>
          </div>
          {/* Visualisations */}
          <Tabs value={activeTab} onValueChange={setActiveTab}>
            <TabsList>
              <TabsTrigger value="overview">Overview</TabsTrigger>
              <TabsTrigger value="performers">Top Performers</TabsTrigger>
              <TabsTrigger value="mutations">Mutation Fingerprint</TabsTrigger>
              <TabsTrigger value="landscape">Activity Landscape</TabsTrigger>
            </TabsList>

            <TabsContent value="overview" className="space-y-6 mt-6">
              <ActivityDistributionChart
                variants={passedVariants}
                experimentId={selectedExpId}
              />

              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center gap-2">
                    <TrendingUp className="h-5 w-5" />
                    Top 10 Performing Variants
                  </CardTitle>
                  <CardDescription>
                    Variants with highest activity scores
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <TopPerformersTable
                    variants={topPerformers}
                    onSelectVariant={(idx) => {
                      setSelectedVariantIndex(idx);
                      setActiveTab("mutations");
                    }}
                    selectedVariant={selectedVariantIndex}
                    isLoading={isLoadingTopPerformers}
                  />
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="performers" className="mt-6">
              <Card>
                <CardHeader>
                  <CardTitle>Top Performing Variants</CardTitle>
                  <CardDescription>
                    Detailed view of the top 10 variants by unified activity
                    score
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <TopPerformersTable
                    variants={topPerformers}
                    onSelectVariant={(idx) => {
                      setSelectedVariantIndex(idx);
                      setActiveTab("mutations");
                    }}
                    selectedVariant={selectedVariantIndex}
                    isLoading={isLoadingTopPerformers}
                    detailed
                  />
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="mutations" className="mt-6">
              <MutationFingerprint
                variants={topPerformers}
                experimentId={selectedExpId}
                selectedVariantIndex={selectedVariantIndex}
                onSelectVariant={(idx) => setSelectedVariantIndex(idx)}
              />
            </TabsContent>

            <TabsContent value="landscape" className="mt-6">
              <ActivityLandscape
                variants={passedVariants}
                experimentId={selectedExpId}
              />
            </TabsContent>
          </Tabs>
        </>
      )}
    </div>
  );
}

export default function AnalysisPage() {
  return (
    <Suspense fallback={<Loading />}>
      <AnalysisContent />
    </Suspense>
  );
}
