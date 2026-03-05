"use client";

import React from "react";

import { useEffect, useState, useRef } from "react";
import Link from "next/link";
import { useRouter, useParams } from "next/navigation";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from "@/components/ui/select";
import { toast } from "@/hooks/use-toast";
import {
  ArrowLeft,
  Upload,
  FileText,
  CheckCircle2,
  XCircle,
  Loader2,
  BarChart3,
  Dna,
  AlertTriangle,
  X,
  RefreshCw,
} from "lucide-react";

import type { Experiment, VariantData } from "@/lib/types";

type MappingEntry = { rawCol: string; canonicalField: string };

const FIELD_LABELS: Record<string, string> = {
  plasmid_variant_index: "Plasmid Variant Index",
  parent_plasmid_variant: "Parent Plasmid Variant",
  generation: "Generation",
  assembled_dna_sequence: "Assembled DNA Sequence",
  dna_yield: "DNA Yield",
  protein_yield: "Protein Yield",
  is_control: "Is Control",
};

export default function ExperimentDetailPage() {
  const params = useParams();
  const id = params.id as string;
  const router = useRouter();

  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [variants, setVariants] = useState<VariantData[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [isUploading, setIsUploading] = useState(false);
  const [isAnalyzing, setIsAnalyzing] = useState(false);
  const [justCompleted, setJustCompleted] = useState(false);
  const [analysisBanner, setAnalysisBanner] = useState<{
    message: string;
  } | null>(null);
  // Elapsed-seconds counter shown while analysis is running
  const [elapsedSecs, setElapsedSecs] = useState(0);
  const analysisStartRef = useRef<number | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);
  const [uploadResult, setUploadResult] = useState<{
    parsed: number;
    processed: number;
    passedQC: number;
    failedQC: number;
    errors: Array<{ row: number; field: string; message: string }>;
    warnings: string[];
  } | null>(null);
  const [uploadError, setUploadError] = useState<string | null>(null);

  // Column-mapping confirmation state
  const [isPreviewingMapping, setIsPreviewingMapping] = useState(false);
  const [showMappingReview, setShowMappingReview] = useState(false);
  const [pendingContent, setPendingContent] = useState("");
  const [pendingFormat, setPendingFormat] = useState("tsv");
  const [mappingEntries, setMappingEntries] = useState<MappingEntry[]>([]);
  const [canonicalFields, setCanonicalFields] = useState<string[]>([]);

  useEffect(() => {
    loadExperiment(); // fast — clears the full-page spinner
    loadVariants(); // slow — populates variant-dependent UI in background
  }, [id]);

  // Tick the elapsed-seconds counter while analysis is running
  useEffect(() => {
    if (experiment?.analysisStatus !== "analyzing") {
      analysisStartRef.current = null;
      setElapsedSecs(0);
      return;
    }
    // Record start time on first render showing "analyzing"
    if (analysisStartRef.current === null) {
      analysisStartRef.current = Date.now();
    }
    const tick = setInterval(() => {
      setElapsedSecs(
        Math.floor(
          (Date.now() - (analysisStartRef.current ?? Date.now())) / 1000,
        ),
      );
    }, 1000);
    return () => clearInterval(tick);
  }, [experiment?.analysisStatus]);

  // Manually check whether analysis has finished — called from the banner button
  const [isCheckingStatus, setIsCheckingStatus] = useState(false);
  const checkAnalysisStatus = async () => {
    setIsCheckingStatus(true);
    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${id}?include_variants=false`,
        { credentials: "include" },
      );
      const data = await res.json();
      if (data.success) {
        setExperiment(data.experiment);
        if (data.experiment.analysisStatus === "completed") {
          setIsAnalyzing(false);
          setJustCompleted(true);
          setTimeout(() => setJustCompleted(false), 10000);
          setAnalysisBanner({
            message:
              data.experiment.analysisMessage ||
              "Sequences analysed successfully",
          });
          loadVariants();
          toast({
            title: "Analysis Complete!",
            description:
              data.experiment.analysisMessage ||
              "Sequences analysed successfully",
          });
        } else if (data.experiment.analysisStatus === "failed") {
          setIsAnalyzing(false);
          toast({
            title: "Analysis Failed",
            description: data.experiment.analysisMessage || "An error occurred",
            variant: "destructive",
          });
        } else {
          toast({
            title: "Still running",
            description: "Analysis is still in progress — check again shortly.",
          });
        }
      }
    } catch {
      toast({ title: "Could not reach server", variant: "destructive" });
    } finally {
      setIsCheckingStatus(false);
    }
  };

  // Fast: metadata only — gates the initial page display
  const loadExperiment = async () => {
    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${id}?include_variants=false`,
        { credentials: "include" },
      );
      const data = await res.json();
      if (data.success) {
        setExperiment(data.experiment);
      } else {
        toast({ title: "Experiment not found", variant: "destructive" });
        router.push("/dashboard/experiments");
      }
    } catch (error) {
      console.error("Failed to load experiment:", error);
      toast({ title: "Failed to load experiment", variant: "destructive" });
    } finally {
      setIsLoading(false);
    }
  };

  // Slower: loads variant rows in the background — doesn't block page display
  const loadVariants = async () => {
    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${id}`,
        { credentials: "include" },
      );
      const data = await res.json();
      if (data.success) {
        setVariants(data.variants || []);
      }
    } catch {
      // Silently ignore — not critical for page structure
    }
  };

  // Step 1: read file + fetch auto-detected mapping from backend
  const handleFileChange = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setIsPreviewingMapping(true);
    setUploadResult(null);
    setUploadError(null);
    setShowMappingReview(false);

    try {
      const content = await file.text();
      const format = file.name.endsWith(".json") ? "json" : "tsv";

      setPendingContent(content);
      setPendingFormat(format);

      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${id}/preview-mapping`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          credentials: "include",
          body: JSON.stringify({ data: content, format }),
        },
      );

      const result = await res.json();

      if (result.success) {
        const entries: MappingEntry[] = (result.raw_columns as string[]).map(
          (raw) => ({
            rawCol: raw,
            canonicalField:
              (result.column_mapping as Record<string, string>)[raw] ??
              "__metadata__",
          }),
        );
        setMappingEntries(entries);
        setCanonicalFields((result.canonical_fields as string[]) ?? []);
        setShowMappingReview(true);
      } else {
        toast({
          title: result.error || "Failed to read file",
          variant: "destructive",
        });
      }
    } catch (error) {
      console.error("Preview mapping error:", error);
      toast({ title: "Failed to read file", variant: "destructive" });
    } finally {
      setIsPreviewingMapping(false);
    }
  };

  // Step 2: user confirmed the mapping — send the file for full processing
  const handleConfirmMapping = async () => {
    if (!pendingContent) return;

    setIsUploading(true);
    setShowMappingReview(false);

    const columnMapping = Object.fromEntries(
      mappingEntries
        .filter(
          (e) => e.canonicalField !== "" && e.canonicalField !== "__metadata__",
        )
        .map((e) => [e.rawCol, e.canonicalField]),
    );

    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${id}/upload-data`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          credentials: "include",
          body: JSON.stringify({
            data: pendingContent,
            format: pendingFormat,
            column_mapping: columnMapping,
          }),
        },
      );

      const result = await res.json();

      if (result.success) {
        setUploadResult(result);
        setUploadError(null);
        toast({
          title: "Success!",
          description: `Parsed ${result.processed} variants successfully. ${result.failedQC} rows failed QC.`,
        });
        loadExperiment();
        loadVariants();
      } else {
        const msg = result.error || "Failed to parse data";
        setUploadError(msg);
        toast({ title: msg, variant: "destructive" });
      }
    } catch (error) {
      console.error("Upload error:", error);
      const msg =
        "Failed to upload file — check the browser console for details.";
      setUploadError(msg);
      toast({ title: "Failed to upload file", variant: "destructive" });
    } finally {
      setIsUploading(false);
      setPendingContent("");
      setPendingFormat("tsv");
    }
  };

  const handleCancelMapping = () => {
    setShowMappingReview(false);
    setPendingContent("");
    setPendingFormat("tsv");
    setMappingEntries([]);
    setUploadError(null);
    if (fileInputRef.current) fileInputRef.current.value = "";
  };

  const handleAnalyzeSequences = async () => {
    if (!experiment) return;

    setIsAnalyzing(true);
    setJustCompleted(false);
    setAnalysisBanner(null);

    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments/${id}/analyze-sequences`,
        {
          method: "POST",
          credentials: "include",
        },
      );

      const result = await res.json();

      if (result.success) {
        toast({
          title: "Analysis Started",
          description:
            "Running in background — status will update automatically.",
        });
        setExperiment((prev) =>
          prev ? { ...prev, analysisStatus: "analyzing" } : prev,
        );
        // Clear the POST-in-flight flag immediately.
        // From here, the button loading state is driven purely by
        // experiment.analysisStatus (kept current by the polling effect).
        setIsAnalyzing(false);
      } else {
        toast({
          title: result.error || "Analysis failed",
          variant: "destructive",
        });
        setIsAnalyzing(false);
      }
    } catch (error) {
      console.error("Analysis error:", error);
      toast({ title: "Failed to analyse sequences", variant: "destructive" });
      setIsAnalyzing(false);
    }
  };

  // Derived: which required canonical fields are not yet assigned in the mapping review
  const REQUIRED_CANONICAL_FIELDS = [
    "plasmid_variant_index",
    "generation",
    "assembled_dna_sequence",
    "dna_yield",
    "protein_yield",
    "is_control",
  ];
  const assignedCanonicals = new Set(
    mappingEntries
      .filter((e) => e.canonicalField !== "__metadata__")
      .map((e) => e.canonicalField),
  );
  const missingRequiredFields = REQUIRED_CANONICAL_FIELDS.filter(
    (f) => !assignedCanonicals.has(f),
  );

  if (isLoading) {
    return (
      <div className="flex items-center justify-center py-12">
        <Loader2 className="h-6 w-6 animate-spin text-muted-foreground" />
      </div>
    );
  }

  if (!experiment) {
    return null;
  }

  const passedVariants = variants.filter((v) => v.qcStatus === "passed");
  const generations = [...new Set(variants.map((v) => v.generation))].sort(
    (a, b) => a - b,
  );

  return (
    <div className="space-y-6">
      <div className="flex items-center gap-4">
        <Button variant="ghost" size="icon" asChild>
          <Link href="/dashboard/experiments">
            <ArrowLeft className="h-4 w-4" />
          </Link>
        </Button>
        <div className="flex-1">
          <h1 className="text-2xl font-bold text-foreground">
            {experiment.name}
          </h1>
          <p className="text-muted-foreground">
            {experiment.proteinName || experiment.proteinAccession} (
            {experiment.proteinAccession})
          </p>
        </div>
        <Badge
          variant={
            experiment.validationStatus === "valid"
              ? "default"
              : experiment.validationStatus === "invalid"
                ? "destructive"
                : "secondary"
          }
        >
          {experiment.validationStatus}
        </Badge>
      </div>

      {/* Global analysis status banners — always visible regardless of active tab */}
      {experiment.analysisStatus === "analyzing" && (
        <div className="flex items-center gap-3 p-3 rounded-lg bg-blue-50 border border-blue-200 text-blue-800 dark:bg-blue-950 dark:border-blue-800 dark:text-blue-300">
          <Loader2 className="h-4 w-4 animate-spin shrink-0 text-blue-600 dark:text-blue-400" />
          <p className="text-sm flex-1">
            Sequence analysis in progress
            {elapsedSecs > 0 && (
              <span className="ml-1 font-mono text-xs">
                ({elapsedSecs}s elapsed)
              </span>
            )}
            {" — this may take a few minutes for large datasets."}
          </p>
          <button
            onClick={checkAnalysisStatus}
            disabled={isCheckingStatus}
            className="shrink-0 flex items-center gap-1.5 text-xs font-medium px-2.5 py-1 rounded-md bg-blue-100 hover:bg-blue-200 disabled:opacity-50 text-blue-800 dark:bg-blue-900 dark:hover:bg-blue-800 dark:text-blue-200 transition-colors"
          >
            {isCheckingStatus ? (
              <Loader2 className="h-3 w-3 animate-spin" />
            ) : (
              <RefreshCw className="h-3 w-3" />
            )}
            Check status
          </button>
        </div>
      )}
      {analysisBanner && (
        <div className="flex items-start gap-3 p-4 rounded-lg bg-green-50 border border-green-200 text-green-800 dark:bg-green-950 dark:border-green-800 dark:text-green-300">
          <CheckCircle2 className="h-5 w-5 mt-0.5 shrink-0 text-green-600" />
          <div className="flex-1">
            <p className="text-sm font-semibold">Analysis Complete!</p>
            <p className="text-sm">{analysisBanner.message}</p>
          </div>
          <button
            onClick={() => setAnalysisBanner(null)}
            className="shrink-0 text-green-600 hover:text-green-800 dark:text-green-400 dark:hover:text-green-200"
          >
            <X className="h-4 w-4" />
          </button>
        </div>
      )}
      <Tabs defaultValue="overview">
        <TabsList>
          <TabsTrigger value="overview">Overview</TabsTrigger>
          <TabsTrigger value="upload">Upload Data</TabsTrigger>
          {variants.length > 0 && (
            <TabsTrigger value="variants">
              Variants ({variants.length})
            </TabsTrigger>
          )}
        </TabsList>

        <TabsContent value="overview" className="space-y-6 mt-6">
          {/* Experiment info */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Dna className="h-5 w-5" />
                  Protein Information
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-3">
                <div>
                  <p className="text-sm text-muted-foreground">Name</p>
                  <p className="font-medium">
                    {experiment.proteinName || "Unknown"}
                  </p>
                </div>
                <div>
                  <p className="text-sm text-muted-foreground">Organism</p>
                  <p className="font-medium">
                    {(experiment.proteinFeatures as any)?.organism || "Unknown"}
                  </p>
                </div>
                <div>
                  <p className="text-sm text-muted-foreground">Length</p>
                  <p className="font-mono">
                    {experiment.wtProteinSequence?.length || 0} aa
                  </p>
                </div>
                <div>
                  <p className="text-sm text-muted-foreground">Accession</p>
                  <a
                    href={`https://www.uniprot.org/uniprotkb/${experiment.proteinAccession}`}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="font-mono text-primary hover:underline"
                  >
                    {experiment.proteinAccession}
                  </a>
                </div>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <FileText className="h-5 w-5" />
                  Plasmid Information
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-3">
                <div>
                  <p className="text-sm text-muted-foreground">Name</p>
                  <p className="font-medium">{experiment.plasmidName}</p>
                </div>
                <div>
                  <p className="text-sm text-muted-foreground">Length</p>
                  <p className="font-mono">
                    {
                      experiment.plasmidSequence
                        .split("\n")
                        .filter((line) => !line.startsWith(">"))
                        .join("")
                        .replace(/\s/g, "").length
                    }{" "}
                    bp
                  </p>
                </div>
                <div>
                  <p className="text-sm text-muted-foreground">Validation</p>
                  <div className="flex items-center gap-2 mt-1">
                    {experiment.validationStatus === "valid" ? (
                      <CheckCircle2 className="h-4 w-4 text-accent" />
                    ) : (
                      <XCircle className="h-4 w-4 text-destructive" />
                    )}
                    <p className="text-sm">{experiment.validationMessage}</p>
                  </div>
                </div>
              </CardContent>
            </Card>
          </div>

          {/* Stats if data uploaded */}
          {variants.length > 0 && (
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <Card>
                <CardContent className="pt-6">
                  <p className="text-sm text-muted-foreground">
                    Total Variants
                  </p>
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
                  <p className="text-sm text-muted-foreground">Generations</p>
                  <p className="text-2xl font-bold">{generations.length}</p>
                </CardContent>
              </Card>
              <Card>
                <CardContent className="pt-6">
                  <p className="text-sm text-muted-foreground">
                    Highest Activity
                  </p>
                  <p className="text-2xl font-bold">
                    {Math.max(...variants.map((v) => v.activityScore)).toFixed(
                      2,
                    )}
                  </p>
                </CardContent>
              </Card>
            </div>
          )}

          {variants.length > 0 && (
            <Card>
              <CardContent className="pt-6">
                <Button asChild>
                  <Link href={`/dashboard/analysis?experiment=${id}`}>
                    <BarChart3 className="h-4 w-4 mr-2" />
                    View Full Analysis
                  </Link>
                </Button>
              </CardContent>
            </Card>
          )}
        </TabsContent>

        <TabsContent value="upload" className="mt-6">
          {/* Flash-green keyframe for mapping rows */}
          <style>{`
            @keyframes mapping-flash {
              0%   { background-color: transparent; }
              25%  { background-color: hsl(142 70% 45% / 0.14); }
              100% { background-color: transparent; }
            }
            .mapping-row-flash { animation: mapping-flash 1.1s ease-out forwards; }
          `}</style>

          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Upload className="h-5 w-5" />
                Upload Experimental Data
              </CardTitle>
              <CardDescription>
                {showMappingReview
                  ? "Review and adjust the detected column mapping before uploading."
                  : "Upload your directed evolution data in TSV or JSON format"}
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              {/* ── File picker (hidden while reviewing) ─────────────────── */}
              {!showMappingReview && (
                <div className="border-2 border-dashed border-border rounded-lg p-8 text-center hover:border-primary/50 transition-colors">
                  <input
                    ref={fileInputRef}
                    type="file"
                    accept=".tsv,.json,.txt"
                    onChange={handleFileChange}
                    className="hidden"
                    id="data-upload"
                    disabled={isUploading || isPreviewingMapping}
                  />
                  <label htmlFor="data-upload" className="cursor-pointer">
                    {isUploading || isPreviewingMapping ? (
                      <Loader2 className="h-8 w-8 mx-auto text-muted-foreground mb-2 animate-spin" />
                    ) : (
                      <Upload className="h-8 w-8 mx-auto text-muted-foreground mb-2" />
                    )}
                    <p className="text-sm text-muted-foreground">
                      {isUploading
                        ? "Uploading..."
                        : isPreviewingMapping
                          ? "Reading file..."
                          : "Click to upload TSV or JSON file"}
                    </p>
                    <p className="text-xs text-muted-foreground mt-1">
                      Required fields: Plasmid_Variant_Index,
                      Assembled_DNA_Sequence, DNA_Yield, Protein_Yield
                    </p>
                  </label>
                </div>
              )}

              {/* ── Persistent upload error ──────────────────────────────── */}
              {uploadError && !showMappingReview && (
                <div className="rounded-lg border border-destructive/40 bg-destructive/5 p-4 flex items-start gap-3">
                  <XCircle className="h-5 w-5 text-destructive mt-0.5 shrink-0" />
                  <div className="flex-1 min-w-0">
                    <p className="font-semibold text-destructive text-sm mb-0.5">
                      Upload failed
                    </p>
                    <p className="text-sm text-destructive/90 break-words">
                      {uploadError}
                    </p>
                    <p className="text-xs text-destructive/70 mt-1">
                      Fix your file or adjust the column mapping and try again.
                    </p>
                  </div>
                  <button
                    onClick={() => setUploadError(null)}
                    className="shrink-0 text-destructive/60 hover:text-destructive"
                  >
                    <X className="h-4 w-4" />
                  </button>
                </div>
              )}

              {/* ── Column mapping review ─────────────────────────────────── */}
              {showMappingReview && (
                <div className="space-y-4">
                  <p className="text-sm text-muted-foreground">
                    The system auto-detected the mapping below. Each row flashes
                    green to confirm detection. Adjust any dropdown before
                    confirming.
                  </p>

                  <div className="rounded-lg border border-border overflow-hidden">
                    <table className="w-full text-sm">
                      <thead>
                        <tr className="bg-muted/50 border-b">
                          <th className="text-left px-4 py-2.5 font-medium text-muted-foreground w-1/2">
                            File Column
                          </th>
                          <th className="text-left px-4 py-2.5 font-medium text-muted-foreground w-1/2">
                            Maps To
                          </th>
                        </tr>
                      </thead>
                      <tbody>
                        {mappingEntries.map((entry, idx) => (
                          <tr
                            key={entry.rawCol}
                            className="mapping-row-flash border-b last:border-0"
                            style={{ animationDelay: `${idx * 70}ms` }}
                          >
                            <td className="px-4 py-2.5 font-mono text-xs text-foreground">
                              {entry.rawCol}
                            </td>
                            <td className="px-4 py-2.5">
                              <Select
                                value={entry.canonicalField}
                                onValueChange={(val) =>
                                  setMappingEntries((prev) =>
                                    prev.map((e, i) =>
                                      i === idx
                                        ? { ...e, canonicalField: val }
                                        : e,
                                    ),
                                  )
                                }
                              >
                                <SelectTrigger className="h-8 text-xs">
                                  <SelectValue placeholder="— extra metadata —" />
                                </SelectTrigger>
                                <SelectContent>
                                  <SelectItem value="__metadata__">
                                    — extra metadata —
                                  </SelectItem>
                                  {canonicalFields.map((field) => (
                                    <SelectItem key={field} value={field}>
                                      {FIELD_LABELS[field] ?? field}
                                    </SelectItem>
                                  ))}
                                </SelectContent>
                              </Select>
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>

                  {/* Error if required fields are missing — blocks upload */}
                  {missingRequiredFields.length > 0 && (
                    <div className="rounded-lg border border-destructive/40 bg-destructive/5 p-3 flex items-start gap-2">
                      <XCircle className="h-4 w-4 text-destructive mt-0.5 shrink-0" />
                      <p className="text-sm text-destructive">
                        <span className="font-semibold">
                          Upload blocked — missing required columns:{" "}
                        </span>
                        <span className="font-medium">
                          {missingRequiredFields
                            .map((f) => FIELD_LABELS[f] ?? f)
                            .join(", ")}
                        </span>
                        . Your file must contain these columns (or map them
                        above) before you can upload.
                      </p>
                    </div>
                  )}

                  <div className="flex gap-3">
                    <Button
                      onClick={handleConfirmMapping}
                      disabled={isUploading || missingRequiredFields.length > 0}
                      title={
                        missingRequiredFields.length > 0
                          ? `Upload blocked: missing ${missingRequiredFields.map((f) => FIELD_LABELS[f] ?? f).join(", ")}`
                          : undefined
                      }
                    >
                      {isUploading ? (
                        <>
                          <Loader2 className="h-4 w-4 animate-spin mr-2" />
                          Uploading...
                        </>
                      ) : (
                        <>
                          <CheckCircle2 className="h-4 w-4 mr-2" />
                          Confirm Mapping &amp; Upload
                        </>
                      )}
                    </Button>
                    <Button
                      variant="outline"
                      onClick={handleCancelMapping}
                      disabled={isUploading}
                    >
                      Start Over
                    </Button>
                  </div>
                </div>
              )}

              {uploadResult && (
                <div className="space-y-4">
                  <div className="rounded-lg border border-accent/30 bg-accent/5 p-4">
                    <div className="flex items-center gap-2 mb-3">
                      <CheckCircle2 className="h-5 w-5 text-accent" />
                      <p className="font-medium">Data Processed Successfully</p>
                    </div>
                    <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                      <div>
                        <p className="text-muted-foreground">Rows Parsed</p>
                        <p className="font-medium">{uploadResult.parsed}</p>
                      </div>
                      <div>
                        <p className="text-muted-foreground">
                          Variants Processed
                        </p>
                        <p className="font-medium">{uploadResult.processed}</p>
                      </div>
                      <div>
                        <p className="text-muted-foreground">Passed QC</p>
                        <p className="font-medium text-accent">
                          {uploadResult.passedQC}
                        </p>
                      </div>
                      <div>
                        <p className="text-muted-foreground">Failed QC</p>
                        <p className="font-medium text-destructive">
                          {uploadResult.failedQC}
                        </p>
                      </div>
                    </div>
                  </div>

                  {uploadResult.warnings.length > 0 && (
                    <div className="rounded-lg border border-yellow-500/30 bg-yellow-500/5 p-4">
                      <div className="flex items-center gap-2 mb-2">
                        <AlertTriangle className="h-4 w-4 text-yellow-600" />
                        <p className="font-medium text-yellow-600">Warnings</p>
                      </div>
                      <ul className="text-sm text-yellow-600/80 space-y-1">
                        {uploadResult.warnings.map((w, i) => (
                          <li key={i}>{w}</li>
                        ))}
                      </ul>
                    </div>
                  )}

                  {uploadResult.errors.length > 0 && (
                    <div className="rounded-lg border border-destructive/30 bg-destructive/5 p-4">
                      <div className="flex items-center gap-2 mb-2">
                        <XCircle className="h-4 w-4 text-destructive" />
                        <p className="font-medium text-destructive">
                          QC Errors ({uploadResult.errors.length})
                        </p>
                      </div>
                      <ul className="text-sm text-destructive/80 space-y-1 max-h-40 overflow-y-auto">
                        {uploadResult.errors.slice(0, 10).map((e, i) => (
                          <li key={i}>
                            Row {e.row}: {e.field} - {e.message}
                          </li>
                        ))}
                        {uploadResult.errors.length > 10 && (
                          <li>
                            ... and {uploadResult.errors.length - 10} more
                          </li>
                        )}
                      </ul>
                    </div>
                  )}
                </div>
              )}
            </CardContent>
          </Card>

          {/* Sequence Analysis Section */}
          {variants.length > 0 && (
            <Card className="mt-6">
              <CardHeader>
                <CardTitle className="flex items-center gap-2">
                  <Dna className="h-5 w-5" />
                  Sequence Analysis
                </CardTitle>
                <CardDescription>
                  Analyse DNA sequences to detect protein mutations
                </CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="flex items-center justify-between p-4 border rounded-lg">
                  <div className="flex-1">
                    <p className="font-medium mb-1">Analysis Status</p>
                    <div className="flex items-center gap-2">
                      {experiment.analysisStatus === "completed" ? (
                        <>
                          <CheckCircle2 className="h-4 w-4 text-accent" />
                          <span className="text-sm text-muted-foreground">
                            {experiment.analysisMessage || "Analysis completed"}
                          </span>
                        </>
                      ) : experiment.analysisStatus === "analyzing" ? (
                        <>
                          <Loader2 className="h-4 w-4 animate-spin text-primary" />
                          <span className="text-sm text-muted-foreground">
                            Analysis in progress — see banner above
                          </span>
                        </>
                      ) : experiment.analysisStatus === "failed" ? (
                        <>
                          <XCircle className="h-4 w-4 text-destructive" />
                          <span className="text-sm text-destructive">
                            {experiment.analysisMessage || "Analysis failed"}
                          </span>
                        </>
                      ) : (
                        <>
                          <AlertTriangle className="h-4 w-4 text-yellow-500" />
                          <span className="text-sm text-muted-foreground">
                            Not yet analyzed
                          </span>
                        </>
                      )}
                    </div>
                  </div>
                  <Button
                    onClick={handleAnalyzeSequences}
                    disabled={
                      isAnalyzing || experiment.analysisStatus === "analyzing"
                    }
                    variant="outline"
                    className={
                      justCompleted
                        ? "border-green-500 bg-green-50 text-green-700 hover:bg-green-100 dark:bg-green-950 dark:text-green-400 dark:border-green-700 dark:hover:bg-green-900"
                        : experiment.analysisStatus === "completed"
                          ? "border-green-500 text-green-700 hover:bg-green-50 dark:text-green-400 dark:border-green-700 dark:hover:bg-green-950"
                          : "border-transparent bg-primary text-primary-foreground hover:bg-primary/90"
                    }
                  >
                    {isAnalyzing ||
                    experiment.analysisStatus === "analyzing" ? (
                      <>
                        <Loader2 className="h-4 w-4 mr-2 animate-spin" />
                        Analysing…
                      </>
                    ) : justCompleted ? (
                      <>
                        <CheckCircle2 className="h-4 w-4 mr-2" />
                        Analysis Complete!
                      </>
                    ) : experiment.analysisStatus === "completed" ? (
                      <>
                        <CheckCircle2 className="h-4 w-4 mr-2" />
                        Re-analyse
                      </>
                    ) : (
                      <>
                        <Dna className="h-4 w-4 mr-2" />
                        Analyse Sequences
                      </>
                    )}
                  </Button>
                </div>
                {/* analysisBanner is now shown globally above the tabs */}
                <div className="text-sm text-muted-foreground space-y-1">
                  <p>• Translates DNA sequences to protein sequences</p>
                  <p>• Detects mutations compared to wild-type</p>
                  <p>• Analysis may take several minutes for large datasets</p>
                </div>
              </CardContent>
            </Card>
          )}
        </TabsContent>

        {variants.length > 0 && (
          <TabsContent value="variants" className="mt-6">
            <Card>
              <CardHeader>
                <CardTitle>Variant Data</CardTitle>
                <CardDescription>
                  Showing first 50 variants sorted by activity score
                </CardDescription>
              </CardHeader>
              <CardContent>
                <div className="overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b">
                        <th className="text-left py-2 px-2">Index</th>
                        <th className="text-left py-2 px-2">Gen</th>
                        <th className="text-left py-2 px-2">DNA Yield</th>
                        <th className="text-left py-2 px-2">Protein</th>
                        <th className="text-left py-2 px-2">Activity</th>
                        <th className="text-left py-2 px-2">QC</th>
                      </tr>
                    </thead>
                    <tbody>
                      {variants
                        .sort((a, b) => b.activityScore - a.activityScore)
                        .slice(0, 50)
                        .map((v) => (
                          <tr key={v.id} className="border-b">
                            <td className="py-2 px-2 font-mono">
                              {v.plasmidVariantIndex}
                            </td>
                            <td className="py-2 px-2">{v.generation}</td>
                            <td className="py-2 px-2">
                              {v.dnaYield.toFixed(2)}
                            </td>
                            <td className="py-2 px-2">
                              {v.proteinYield.toFixed(2)}
                            </td>
                            <td className="py-2 px-2 font-medium">
                              {v.activityScore.toFixed(3)}
                            </td>
                            <td className="py-2 px-2">
                              <Badge
                                variant={
                                  v.qcStatus === "passed"
                                    ? "default"
                                    : "destructive"
                                }
                              >
                                {v.qcStatus}
                              </Badge>
                            </td>
                          </tr>
                        ))}
                    </tbody>
                  </table>
                </div>
              </CardContent>
            </Card>
          </TabsContent>
        )}
      </Tabs>
    </div>
  );
}
