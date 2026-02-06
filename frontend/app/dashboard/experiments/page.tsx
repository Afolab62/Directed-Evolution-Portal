"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
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
import { Badge } from "@/components/ui/badge";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import {
  AlertDialog,
  AlertDialogAction,
  AlertDialogCancel,
  AlertDialogContent,
  AlertDialogDescription,
  AlertDialogFooter,
  AlertDialogHeader,
  AlertDialogTitle,
} from "@/components/ui/alert-dialog";
import { FlaskConical, Plus, ChevronRight, Trash2 } from "lucide-react";
import { toast } from "sonner";
import type { Experiment } from "@/lib/types";

export default function ExperimentsPage() {
  const { user } = useAuth();
  const router = useRouter();
  const [experiments, setExperiments] = useState<Experiment[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [experimentToDelete, setExperimentToDelete] = useState<{
    id: string;
    name: string;
  } | null>(null);

  const loadExperiments = async () => {
    if (!user) return;

    try {
      const res = await fetch(`/api/experiments?userId=${user.id}`);
      const data = await res.json();
      if (data.success) {
        setExperiments(data.experiments);
      }
    } catch (err) {
      console.error("Failed to load experiments:", err);
    } finally {
      setIsLoading(false);
    }
  };

  useEffect(() => {
    loadExperiments();
  }, [user]);

  const handleDelete = async (id: string, name: string) => {
    setExperimentToDelete({ id, name });
    setDeleteDialogOpen(true);
  };

  const confirmDelete = async () => {
    if (!experimentToDelete) return;

    try {
      const res = await fetch(`/api/experiments/${experimentToDelete.id}`, {
        method: "DELETE",
      });
      const data = await res.json();

      if (data.success) {
        toast.success("Experiment deleted");
        loadExperiments();
      } else {
        toast.error("Failed to delete experiment");
      }
    } catch {
      toast.error("Failed to delete experiment");
    } finally {
      setDeleteDialogOpen(false);
      setExperimentToDelete(null);
    }
  };

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between gap-4">
        <div>
          <h1 className="text-2xl font-bold text-foreground">Experiments</h1>
          <p className="text-muted-foreground mt-1">
            Manage your directed evolution experiments
          </p>
        </div>
        <Button asChild>
          <Link href="/dashboard/new-experiment">
            <Plus className="h-4 w-4 mr-2" />
            New Experiment
          </Link>
        </Button>
      </div>

      <Card>
        <CardHeader>
          <CardTitle>All Experiments</CardTitle>
          <CardDescription>
            Click on an experiment to view details and upload data
          </CardDescription>
        </CardHeader>
        <CardContent>
          {isLoading ? (
            <div className="text-center py-8 text-muted-foreground">
              Loading...
            </div>
          ) : experiments.length === 0 ? (
            <div className="text-center py-12">
              <FlaskConical className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
              <p className="text-muted-foreground mb-4">No experiments yet</p>
              <Button asChild>
                <Link href="/dashboard/new-experiment">
                  <Plus className="h-4 w-4 mr-2" />
                  Create your first experiment
                </Link>
              </Button>
            </div>
          ) : (
            <div className="overflow-x-auto">
              <Table>
                <TableHeader>
                  <TableRow>
                    <TableHead>Name</TableHead>
                    <TableHead>Protein</TableHead>
                    <TableHead>Status</TableHead>
                    <TableHead>Created</TableHead>
                    <TableHead className="text-right">Actions</TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {experiments.map((exp) => (
                    <TableRow
                      key={exp.id}
                      className="cursor-pointer hover:bg-muted/50"
                      onClick={() =>
                        router.push(`/dashboard/experiments/${exp.id}`)
                      }
                    >
                      <TableCell className="font-medium">{exp.name}</TableCell>
                      <TableCell>
                        <div className="text-sm">
                          <p>{exp.proteinName || "Unknown"}</p>
                          <p className="text-muted-foreground font-mono text-xs">
                            {exp.proteinAccession}
                          </p>
                        </div>
                      </TableCell>
                      <TableCell>
                        <Badge
                          variant={
                            exp.validationStatus === "valid"
                              ? "default"
                              : exp.validationStatus === "invalid"
                                ? "destructive"
                                : "secondary"
                          }
                        >
                          {exp.validationStatus}
                        </Badge>
                      </TableCell>
                      <TableCell className="text-muted-foreground text-sm">
                        {new Date(exp.createdAt).toLocaleDateString()}
                      </TableCell>
                      <TableCell className="text-right">
                        <div className="flex items-center justify-end gap-2">
                          <Button
                            variant="ghost"
                            size="icon"
                            onClick={(e) => {
                              e.stopPropagation();
                              handleDelete(exp.id, exp.name);
                            }}
                            className="text-muted-foreground hover:text-destructive"
                          >
                            <Trash2 className="h-4 w-4" />
                            <span className="sr-only">Delete</span>
                          </Button>
                          <ChevronRight className="h-4 w-4 text-muted-foreground" />
                        </div>
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </div>
          )}
        </CardContent>
      </Card>

      <AlertDialog open={deleteDialogOpen} onOpenChange={setDeleteDialogOpen}>
        <AlertDialogContent>
          <AlertDialogHeader>
            <AlertDialogTitle>Delete Experiment</AlertDialogTitle>
            <AlertDialogDescription>
              Are you sure you want to delete &quot;{experimentToDelete?.name}
              &quot;? This action cannot be undone.
            </AlertDialogDescription>
          </AlertDialogHeader>
          <AlertDialogFooter>
            <AlertDialogCancel>Cancel</AlertDialogCancel>
            <AlertDialogAction
              onClick={confirmDelete}
              className="bg-destructive text-destructive-foreground hover:bg-destructive/90"
            >
              Delete
            </AlertDialogAction>
          </AlertDialogFooter>
        </AlertDialogContent>
      </AlertDialog>
    </div>
  );
}
