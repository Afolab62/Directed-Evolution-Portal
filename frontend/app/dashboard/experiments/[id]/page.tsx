"use client";

import React from "react";

import { useEffect, useState } from "react";
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
} from "lucide-react";
import type { Experiment, VariantData } from "@/lib/types";
import { use } from "react"; // Import the use hook

export default function ExperimentDetailPage() {
  const params = useParams();
  const id = params.id as string;
  const router = useRouter();

  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [variants, setVariants] = useState<VariantData[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [isUploading, setIsUploading] = useState(false);
  const [uploadResult, setUploadResult] = useState<{
    parsed: number;
    processed: number;
    passedQC: number;
    failedQC: number;
    errors: Array<{ row: number; field: string; message: string }>;
    warnings: string[];
  } | null>(null);

  useEffect(() => {
    loadExperiment();
  }, [id]);

  const loadExperiment = async () => {
    try {
      const res = await fetch(`${process.env.NEXT_PUBLIC_BACKEND_URL || 'http://localhost:8000'}/api/experiments/${id}`, {
        credentials: 'include'
      });
      const data = await res.json();

      if (data.success) {
        setExperiment(data.experiment);
        setVariants(data.variants || []);
      } else {
        toast({ title: "Experiment not found", variant: "destructive" });
        router.push("/dashboard/experiments");
      }
    } catch {
      toast({ title: "Failed to load experiment", variant: "destructive" });
    } finally {
      setIsLoading(false);
    }
  };

  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    setIsUploading(true);
    setUploadResult(null);

    try {
      // TODO: Implement client-side parsing or backend parse endpoint
      toast({ 
        title: 'Parsing not yet implemented',
        description: 'Parse functionality needs to be migrated',
        variant: 'destructive'
      });
      setIsUploading(false);
      return;
      
      const content = await file.text();
      const format = file.name.endsWith(".json") ? "json" : "tsv";

      const res = await fetch(`/api/experiments/${id}/parse`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ data: content, format }),
      });

      const result = await res.json();

      if (result.success) {
        setUploadResult(result);
        toast({ title: `Parsed ${result.processed} variants successfully` });
        loadExperiment();
      } else {
        toast({
          title: result.error || "Failed to parse data",
          variant: "destructive",
        });
      }
    } catch {
      toast({ title: "Failed to upload file", variant: "destructive" });
    } finally {
      setIsUploading(false);
    }
  };

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
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Upload className="h-5 w-5" />
                Upload Experimental Data
              </CardTitle>
              <CardDescription>
                Upload your directed evolution data in TSV or JSON format
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="border-2 border-dashed border-border rounded-lg p-8 text-center hover:border-primary/50 transition-colors">
                <input
                  type="file"
                  accept=".tsv,.json,.txt"
                  onChange={handleFileUpload}
                  className="hidden"
                  id="data-upload"
                  disabled={isUploading}
                />
                <label htmlFor="data-upload" className="cursor-pointer">
                  {isUploading ? (
                    <Loader2 className="h-8 w-8 mx-auto text-muted-foreground mb-2 animate-spin" />
                  ) : (
                    <Upload className="h-8 w-8 mx-auto text-muted-foreground mb-2" />
                  )}
                  <p className="text-sm text-muted-foreground">
                    {isUploading
                      ? "Processing..."
                      : "Click to upload TSV or JSON file"}
                  </p>
                  <p className="text-xs text-muted-foreground mt-1">
                    Required fields: Plasmid_Variant_Index,
                    Assembled_DNA_Sequence, DNA_Yield, Protein_Yield
                  </p>
                </label>
              </div>

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
                        <th className="text-left py-2 px-2">Mutations</th>
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
                            <td className="py-2 px-2">{v.mutations.length}</td>
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
