"use client";

import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Skeleton } from "@/components/ui/skeleton";
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from "@/components/ui/table";
import { cn } from "@/lib/utils";
import type { VariantData } from "@/lib/types";

interface TopPerformersTableProps {
  variants: VariantData[];
  onSelectVariant?: (index: number) => void;
  selectedVariant?: number | null;
  detailed?: boolean;
  isLoading?: boolean;
}

export function TopPerformersTable({
  variants,
  onSelectVariant,
  selectedVariant,
  detailed = false,
  isLoading = false,
}: TopPerformersTableProps) {
  const columnCount = detailed ? 8 : 6;

  if (isLoading) {
    return (
      <div className="overflow-x-auto">
        <Table>
          <TableHeader>
            <TableRow>
              <TableHead className="w-12">Rank</TableHead>
              <TableHead>Variant Index</TableHead>
              <TableHead>Generation</TableHead>
              <TableHead>Activity Score</TableHead>
              <TableHead>Total Mutations</TableHead>
              {detailed && (
                <>
                  <TableHead>DNA Yield</TableHead>
                  <TableHead>Protein Yield</TableHead>
                </>
              )}
              <TableHead className="text-right">Actions</TableHead>
            </TableRow>
          </TableHeader>
          <TableBody>
            {Array.from({ length: 10 }).map((_, i) => (
              <TableRow key={i}>
                {Array.from({ length: columnCount }).map((_, j) => (
                  <TableCell key={j}>
                    <Skeleton className="h-4 w-full" />
                  </TableCell>
                ))}
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </div>
    );
  }

  if (variants.length === 0) {
    return (
      <div className="text-center py-8 text-muted-foreground">
        No variants to display
      </div>
    );
  }

  return (
    <div className="overflow-x-auto">
      <Table>
        <TableHeader>
          <TableRow>
            <TableHead className="w-12">Rank</TableHead>
            <TableHead>Variant Index</TableHead>
            <TableHead>Generation</TableHead>
            <TableHead>Activity Score</TableHead>
            <TableHead title="Non-synonymous mutations only (amino acid changes)">
              Total Mutations
            </TableHead>
            {detailed && (
              <>
                <TableHead>DNA Yield</TableHead>
                <TableHead>Protein Yield</TableHead>
              </>
            )}
            <TableHead className="text-right">Actions</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {variants.map((variant, idx) => (
            <TableRow
              key={variant.id}
              className={cn(
                selectedVariant === variant.plasmidVariantIndex &&
                  "bg-primary/5",
              )}
            >
              <TableCell className="font-medium">
                {idx === 0 ? (
                  <Badge className="bg-yellow-500 text-yellow-950">#1</Badge>
                ) : idx === 1 ? (
                  <Badge variant="secondary">#2</Badge>
                ) : idx === 2 ? (
                  <Badge variant="outline">#3</Badge>
                ) : (
                  <span className="text-muted-foreground">#{idx + 1}</span>
                )}
              </TableCell>
              <TableCell className="font-mono">
                {variant.plasmidVariantIndex}
              </TableCell>
              <TableCell>
                <Badge variant="outline">Gen {variant.generation}</Badge>
              </TableCell>
              <TableCell className="font-bold text-primary">
                {variant.activityScore.toFixed(3)}
              </TableCell>
              <TableCell>
                <span
                  className={cn(
                    "font-medium",
                    (variant.mutations?.filter((m) => m.type !== "synonymous")
                      .length ?? 0) > 5 && "text-destructive",
                    (variant.mutations?.filter((m) => m.type !== "synonymous")
                      .length ?? 0) > 0 &&
                      (variant.mutations?.filter((m) => m.type !== "synonymous")
                        .length ?? 0) <= 5 &&
                      "text-accent",
                  )}
                >
                  {variant.mutations?.filter((m) => m.type !== "synonymous")
                    .length ?? 0}
                </span>
              </TableCell>
              {detailed && (
                <>
                  <TableCell>{variant.dnaYield.toFixed(2)}</TableCell>
                  <TableCell>{variant.proteinYield.toFixed(2)}</TableCell>
                </>
              )}
              <TableCell className="text-right">
                <Button
                  variant="ghost"
                  size="sm"
                  onClick={() => onSelectVariant?.(variant.plasmidVariantIndex)}
                  className={cn(
                    selectedVariant === variant.plasmidVariantIndex &&
                      "bg-primary/10",
                  )}
                >
                  View Mutations
                </Button>
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </div>
  );
}
