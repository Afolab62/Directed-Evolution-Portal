"use client";

import React, { useState, Suspense } from "react";
import { useRouter } from "next/navigation";
import { useAuth } from "@/lib/auth-context";
import {
  Card,
  CardContent,
  CardDescription,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Badge } from "@/components/ui/badge";
import { toast } from "@/hooks/use-toast";
import {
  Search,
  Upload,
  CheckCircle2,
  XCircle,
  Loader2,
  ArrowRight,
  FileText,
  Dna,
  Download,
} from "lucide-react";
import type { UniProtProtein } from "@/lib/types";
import Loading from "./loading";

type Step = "protein" | "plasmid" | "confirm";

export default function NewExperimentPage() {
  const router = useRouter();
  const { user } = useAuth();

  const [step, setStep] = useState<Step>("protein");
  const [experimentName, setExperimentName] = useState("");
  const [accession, setAccession] = useState("");
  const [protein, setProtein] = useState<UniProtProtein | null>(null);
  const [isSearching, setIsSearching] = useState(false);
  const [searchError, setSearchError] = useState("");

  const [plasmidName, setPlasmidName] = useState("");
  const [plasmidSequence, setPlasmidSequence] = useState("");
  const [isValidating, setIsValidating] = useState(false);
  const [validationResult, setValidationResult] = useState<{
    isValid: boolean;
    message: string;
  } | null>(null);

  const handleSearch = async () => {
    if (!accession.trim()) {
      setSearchError("Please enter a UniProt accession ID");
      return;
    }

    setIsSearching(true);
    setSearchError("");
    setProtein(null);

    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/uniprot/${accession.trim()}`,
        {
          credentials: "include",
        },
      );
      const data = await res.json();

      if (data.success) {
        setProtein(data.protein);
        toast({ title: "Protein found!" });
      } else {
        setSearchError(data.error || "Protein not found");
      }
    } catch {
      setSearchError("Failed to search UniProt");
    } finally {
      setIsSearching(false);
    }
  };

  const handleFileUpload = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      const content = event.target?.result as string;

      // Validate: must have at least one FASTA header line
      const lines = content.split("\n");
      const firstNonEmpty =
        lines.find((l) => l.trim().length > 0)?.trim() ?? "";
      if (!firstNonEmpty.startsWith(">")) {
        toast({
          title: "Invalid file format",
          description:
            'The file does not appear to be a FASTA file. The first line must start with ">" (e.g. ">MyPlasmid"). Please upload a valid .fasta or .fa file.',
          variant: "destructive",
        });
        // Reset the input so the same file can be re-selected after correction
        e.target.value = "";
        return;
      }

      // Validate: sequence lines must only contain valid DNA characters
      const invalidChars = new Set<string>();
      for (const line of lines) {
        if (line.trim().startsWith(">") || line.trim().length === 0) continue;
        for (const ch of line.trim().toUpperCase()) {
          if (!"ACGTN".includes(ch)) invalidChars.add(ch);
        }
      }
      if (invalidChars.size > 0) {
        toast({
          title: "Invalid sequence characters",
          description: `The file contains non-DNA characters: ${[...invalidChars].join(", ")}. DNA sequences may only contain A, C, G, T and N.`,
          variant: "destructive",
        });
        e.target.value = "";
        return;
      }

      setPlasmidSequence(content);

      // Extract name from FASTA header if present
      const headerMatch = content.match(/^>(.+)/m);
      if (headerMatch && !plasmidName) {
        setPlasmidName(headerMatch[1].split(" ")[0]);
      }
    };
    reader.readAsText(file);
  };

  const handleValidate = async () => {
    if (!plasmidSequence.trim()) {
      toast({
        title: "Please upload or paste a plasmid sequence",
        variant: "destructive",
      });
      return;
    }

    setIsValidating(true);
    setValidationResult(null);

    try {
      const res = await fetch(
        `${process.env.NEXT_PUBLIC_BACKEND_URL || "http://localhost:8000"}/api/experiments`,
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          credentials: "include",
          body: JSON.stringify({
            userId: user?.id,
            name: experimentName || `${protein?.name} Experiment`,
            proteinAccession: accession,
            protein,
            plasmidSequence,
            plasmidName,
          }),
        },
      );

      const data = await res.json();

      if (data.success) {
        setValidationResult({
          isValid: data.validation.isValid,
          message: data.validation.message,
        });

        if (data.validation.isValid) {
          toast({ title: "Experiment created successfully!" });
          setStep("confirm");
        } else {
          toast({ title: "Plasmid validation failed", variant: "destructive" });
        }
      } else {
        toast({
          title: data.error || "Failed to create experiment",
          variant: "destructive",
        });
      }
    } catch {
      toast({ title: "Failed to validate plasmid", variant: "destructive" });
    } finally {
      setIsValidating(false);
    }
  };

  return (
    <Suspense fallback={<Loading />}>
      <div className="max-w-3xl mx-auto space-y-6">
        <div>
          <h1 className="text-2xl font-bold text-foreground">New Experiment</h1>
          <p className="text-muted-foreground mt-1">
            Stage a new directed evolution experiment by selecting a protein and
            uploading a plasmid.
          </p>
        </div>

        {/* Progress steps */}
        <div className="flex items-center gap-2">
          <StepIndicator
            number={1}
            label="Select Protein"
            active={step === "protein"}
            completed={step !== "protein"}
          />
          <div className="flex-1 h-px bg-border" />
          <StepIndicator
            number={2}
            label="Upload Plasmid"
            active={step === "plasmid"}
            completed={step === "confirm"}
          />
          <div className="flex-1 h-px bg-border" />
          <StepIndicator
            number={3}
            label="Complete"
            active={step === "confirm"}
            completed={false}
          />
        </div>

        {/* Step 1: Protein Selection */}
        {step === "protein" && (
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Search className="h-5 w-5" />
                Select Wild-Type Protein
              </CardTitle>
              <CardDescription>
                Enter a UniProt accession ID to fetch the protein sequence and
                features. For example: O34996 (BSU DNA Polymerase I)
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="space-y-2">
                <Label htmlFor="experiment-name">
                  Experiment Name (optional)
                </Label>
                <Input
                  id="experiment-name"
                  placeholder="e.g., BSU Pol I Optimization Round 1"
                  value={experimentName}
                  onChange={(e) => setExperimentName(e.target.value)}
                />
              </div>

              <div className="space-y-2">
                <Label htmlFor="accession">UniProt Accession ID</Label>
                <div className="flex gap-2">
                  <Input
                    id="accession"
                    placeholder="e.g., O34996"
                    value={accession}
                    onChange={(e) => setAccession(e.target.value.toUpperCase())}
                    onKeyDown={(e) => e.key === "Enter" && handleSearch()}
                  />
                  <Button onClick={handleSearch} disabled={isSearching}>
                    {isSearching ? (
                      <Loader2 className="h-4 w-4 animate-spin" />
                    ) : (
                      <Search className="h-4 w-4" />
                    )}
                    <span className="ml-2 hidden sm:inline">Search</span>
                  </Button>
                </div>
                {searchError && (
                  <p className="text-sm text-destructive">{searchError}</p>
                )}
              </div>

              {protein && (
                <div className="rounded-lg border border-border p-4 space-y-3 bg-muted/30">
                  <div className="flex items-start justify-between gap-4">
                    <div>
                      <h3 className="font-medium text-foreground">
                        {protein.name}
                      </h3>
                      <p className="text-sm text-muted-foreground">
                        {protein.organism}
                      </p>
                    </div>
                    <Badge variant="secondary">{protein.accession}</Badge>
                  </div>

                  <div className="text-sm">
                    <span className="text-muted-foreground">
                      Sequence length:
                    </span>{" "}
                    <span className="font-mono">{protein.length} aa</span>
                  </div>

                  {protein.features.length > 0 && (
                    <div className="space-y-2">
                      <p className="text-sm font-medium text-muted-foreground">
                        Key Features:
                      </p>
                      <div className="flex flex-wrap gap-1">
                        {protein.features.slice(0, 6).map((feature, idx) => (
                          <Badge
                            key={idx}
                            variant="outline"
                            className="text-xs"
                          >
                            {feature.type}: {feature.location.start}-
                            {feature.location.end}
                          </Badge>
                        ))}
                      </div>
                    </div>
                  )}

                  <div className="pt-2">
                    <Button onClick={() => setStep("plasmid")}>
                      Continue
                      <ArrowRight className="h-4 w-4 ml-2" />
                    </Button>
                  </div>
                </div>
              )}
            </CardContent>
          </Card>
        )}

        {/* Step 2: Plasmid Upload */}
        {step === "plasmid" && (
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Upload className="h-5 w-5" />
                Upload Plasmid Sequence
              </CardTitle>
              <CardDescription>
                Upload the FASTA file containing your plasmid DNA sequence that
                encodes the wild-type protein.
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              {/* Selected protein summary */}
              <div className="rounded-lg border border-border p-3 bg-muted/30 flex items-center gap-3">
                <Dna className="h-5 w-5 text-primary" />
                <div className="flex-1 min-w-0">
                  <p className="font-medium text-sm truncate">
                    {protein?.name}
                  </p>
                  <p className="text-xs text-muted-foreground">
                    {protein?.accession}
                  </p>
                </div>
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={() => setStep("protein")}
                >
                  Change
                </Button>
              </div>

              <div className="space-y-2">
                <Label htmlFor="plasmid-name">Plasmid Name</Label>
                <Input
                  id="plasmid-name"
                  placeholder="e.g., pET28a-BSU_Pol_I_WT"
                  value={plasmidName}
                  onChange={(e) => setPlasmidName(e.target.value)}
                />
              </div>

              <div className="space-y-2">
                <Label>FASTA File Upload</Label>
                <div className="border-2 border-dashed border-border rounded-lg p-6 text-center hover:border-primary/50 transition-colors">
                  <input
                    type="file"
                    accept=".fasta,.fa,.fna,.txt"
                    onChange={handleFileUpload}
                    className="hidden"
                    id="fasta-upload"
                  />
                  <label htmlFor="fasta-upload" className="cursor-pointer">
                    <FileText className="h-8 w-8 mx-auto text-muted-foreground mb-2" />
                    <p className="text-sm text-muted-foreground">
                      Click to upload a FASTA file
                    </p>
                    <p className="text-xs text-muted-foreground mt-1">
                      .fasta, .fa, .fna, or .txt
                    </p>
                  </label>
                </div>
              </div>

              <div className="space-y-2">
                <Label htmlFor="plasmid-sequence">
                  Or paste sequence directly
                </Label>
                <Textarea
                  id="plasmid-sequence"
                  placeholder=">plasmid_name&#10;ATGCGATCG..."
                  value={plasmidSequence}
                  onChange={(e) => setPlasmidSequence(e.target.value)}
                  className="font-mono text-xs min-h-[150px]"
                />
                {plasmidSequence && (
                  <p className="text-xs text-muted-foreground">
                    Sequence length:{" "}
                    {plasmidSequence.replace(/[^ATGCatgc]/g, "").length} bp
                  </p>
                )}
              </div>

              {validationResult && !validationResult.isValid && (
                <div className="rounded-lg border border-destructive/50 bg-destructive/10 p-4">
                  <div className="flex items-start gap-3">
                    <XCircle className="h-5 w-5 text-destructive mt-0.5" />
                    <div>
                      <p className="font-medium text-destructive">
                        Validation Failed
                      </p>
                      <p className="text-sm text-destructive/80 mt-1">
                        {validationResult.message}
                      </p>
                    </div>
                  </div>
                </div>
              )}

              <div className="flex gap-3">
                <Button variant="outline" onClick={() => setStep("protein")}>
                  Back
                </Button>
                <Button
                  onClick={handleValidate}
                  disabled={isValidating || !plasmidSequence.trim()}
                >
                  {isValidating ? (
                    <>
                      <Loader2 className="h-4 w-4 animate-spin mr-2" />
                      Validating...
                    </>
                  ) : (
                    <>
                      Validate & Create
                      <ArrowRight className="h-4 w-4 ml-2" />
                    </>
                  )}
                </Button>
              </div>
            </CardContent>
          </Card>
        )}

        {/* Step 3: Confirmation */}
        {step === "confirm" && (
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2 text-accent">
                <CheckCircle2 className="h-5 w-5" />
                Experiment Created Successfully
              </CardTitle>
              <CardDescription>
                Your experiment has been staged and the plasmid has been
                validated.
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-6">
              <div className="rounded-lg border border-accent/30 bg-accent/5 p-4">
                <p className="text-sm">{validationResult?.message}</p>
              </div>

              <div className="space-y-3">
                <div className="flex justify-between text-sm">
                  <span className="text-muted-foreground">Experiment Name</span>
                  <span className="font-medium">
                    {experimentName || `${protein?.name} Experiment`}
                  </span>
                </div>
                <div className="flex justify-between text-sm">
                  <span className="text-muted-foreground">Protein</span>
                  <span className="font-medium">{protein?.name}</span>
                </div>
                <div className="flex justify-between text-sm">
                  <span className="text-muted-foreground">Accession</span>
                  <span className="font-mono">{protein?.accession}</span>
                </div>
                <div className="flex justify-between text-sm">
                  <span className="text-muted-foreground">Plasmid</span>
                  <span className="font-medium">
                    {plasmidName || "Unnamed"}
                  </span>
                </div>
              </div>

              {/* UniProt Downloads */}
              <div className="rounded-lg border border-border p-4 space-y-3">
                <p className="text-sm font-medium text-foreground flex items-center gap-2">
                  <Download className="h-4 w-4" />
                  Download UniProt Data
                </p>
                <p className="text-xs text-muted-foreground">
                  Download sequence and annotation data for{" "}
                  <span className="font-mono">{protein?.accession}</span> from
                  UniProt.
                </p>
                <div className="flex flex-wrap gap-2">
                  <a
                    href={`${
                      process.env.NEXT_PUBLIC_BACKEND_URL ||
                      "http://localhost:8000"
                    }/api/uniprot/${protein?.accession}/fasta`}
                    download={`${protein?.accession}.fasta`}
                    target="_blank"
                    rel="noreferrer"
                  >
                    <Button variant="outline" size="sm">
                      <Download className="h-3.5 w-3.5 mr-1.5" />
                      FASTA Sequence
                    </Button>
                  </a>
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={() => {
                      if (!protein) return;
                      const blob = new Blob(
                        [JSON.stringify(protein, null, 2)],
                        { type: "application/json" },
                      );
                      const url = URL.createObjectURL(blob);
                      const a = document.createElement("a");
                      a.href = url;
                      a.download = `${protein.accession}_protein_data.json`;
                      a.click();
                      URL.revokeObjectURL(url);
                    }}
                  >
                    <Download className="h-3.5 w-3.5 mr-1.5" />
                    Protein Data (JSON)
                  </Button>
                </div>
              </div>

              <div className="flex gap-3">
                <Button
                  variant="outline"
                  onClick={() => router.push("/dashboard/experiments")}
                >
                  View All Experiments
                </Button>
                <Button
                  onClick={() => router.push("/dashboard/new-experiment")}
                >
                  Create Another
                </Button>
              </div>
            </CardContent>
          </Card>
        )}
      </div>
    </Suspense>
  );
}

function StepIndicator({
  number,
  label,
  active,
  completed,
}: {
  number: number;
  label: string;
  active: boolean;
  completed: boolean;
}) {
  return (
    <div className="flex items-center gap-2">
      <div
        className={`
          w-8 h-8 rounded-full flex items-center justify-center text-sm font-medium
          ${
            completed
              ? "bg-accent text-accent-foreground"
              : active
                ? "bg-primary text-primary-foreground"
                : "bg-muted text-muted-foreground"
          }
        `}
      >
        {completed ? <CheckCircle2 className="h-4 w-4" /> : number}
      </div>
      <span
        className={`text-sm hidden sm:block ${active ? "font-medium text-foreground" : "text-muted-foreground"}`}
      >
        {label}
      </span>
    </div>
  );
}
