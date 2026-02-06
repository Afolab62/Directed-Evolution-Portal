"use client";

import { useEffect, useState } from "react";
import Link from "next/link";
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
  FlaskConical,
  Plus,
  Activity,
  FileCheck,
  AlertCircle,
} from "lucide-react";
import type { Experiment } from "@/lib/types";

export default function DashboardPage() {
  const { user } = useAuth();
  const [experiments, setExperiments] = useState<Experiment[]>([]);
  const [isLoading, setIsLoading] = useState(true);

  useEffect(() => {
    async function loadExperiments() {
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
    }

    loadExperiments();
  }, [user]);

  const validExperiments = experiments.filter(
    (e) => e.validationStatus === "valid",
  );

  return (
    <div className="space-y-8">
      <div>
        <h1 className="text-2xl font-bold text-foreground">Dashboard</h1>
        <p className="text-muted-foreground mt-1">
          Welcome back! Here&apos;s an overview of your directed evolution
          experiments.
        </p>
      </div>

      {/* Stats cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        <Card>
          <CardHeader className="flex flex-row items-center justify-between pb-2">
            <CardTitle className="text-sm font-medium text-muted-foreground">
              Total Experiments
            </CardTitle>
            <FlaskConical className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold">{experiments.length}</div>
          </CardContent>
        </Card>

        <Card>
          <CardHeader className="flex flex-row items-center justify-between pb-2">
            <CardTitle className="text-sm font-medium text-muted-foreground">
              Validated
            </CardTitle>
            <FileCheck className="h-4 w-4 text-muted-foreground" />
          </CardHeader>
          <CardContent>
            <div className="text-2xl font-bold text-accent">
              {validExperiments.length}
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Quick actions */}
      <Card>
        <CardHeader>
          <CardTitle>Quick Actions</CardTitle>
          <CardDescription>
            Get started with your protein engineering workflow
          </CardDescription>
        </CardHeader>
        <CardContent className="flex flex-wrap gap-3">
          <Button asChild>
            <Link href="/dashboard/new-experiment">
              <Plus className="h-4 w-4 mr-2" />
              New Experiment
            </Link>
          </Button>
          <Button variant="outline" asChild>
            <Link href="/dashboard/experiments">
              <FlaskConical className="h-4 w-4 mr-2" />
              View All Experiments
            </Link>
          </Button>
          {validExperiments.length > 0 && (
            <Button variant="outline" asChild>
              <Link href="/dashboard/analysis">
                <Activity className="h-4 w-4 mr-2" />
                View Analysis
              </Link>
            </Button>
          )}
        </CardContent>
      </Card>

      {/* Recent experiments */}
      <Card>
        <CardHeader>
          <CardTitle>Recent Experiments</CardTitle>
          <CardDescription>
            Your latest directed evolution experiments
          </CardDescription>
        </CardHeader>
        <CardContent>
          {isLoading ? (
            <div className="text-center py-8 text-muted-foreground">
              Loading...
            </div>
          ) : experiments.length === 0 ? (
            <div className="text-center py-8">
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
            <div className="space-y-3">
              {experiments.slice(0, 5).map((exp) => (
                <Link
                  key={exp.id}
                  href={`/dashboard/experiments/${exp.id}`}
                  className="block p-4 rounded-lg border border-border hover:bg-muted/50 transition-colors"
                >
                  <div className="flex items-start justify-between gap-4">
                    <div className="min-w-0 flex-1">
                      <h3 className="font-medium text-foreground truncate">
                        {exp.name}
                      </h3>
                      <p className="text-sm text-muted-foreground">
                        {exp.proteinAccession} &middot;{" "}
                        {exp.proteinName || "Unknown protein"}
                      </p>
                    </div>
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
                  </div>
                </Link>
              ))}
            </div>
          )}
        </CardContent>
      </Card>
    </div>
  );
}
