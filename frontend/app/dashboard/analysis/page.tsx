"use client"

import { useEffect, useState, Suspense } from 'react'
import { useSearchParams } from 'next/navigation'
import { useAuth } from '@/lib/auth-context'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Badge } from '@/components/ui/badge'
import { 
  BarChart3, 
  Loader2, 
  TrendingUp, 
  Dna,
  AlertCircle
} from 'lucide-react'
import { TopPerformersTable } from '@/components/analysis/top-performers-table'
import { ActivityDistributionChart } from '@/components/analysis/activity-distribution-chart'
import { MutationFingerprint } from '@/components/analysis/mutation-fingerprint'
import { ActivityLandscape } from '@/components/analysis/activity-landscape'
import type { Experiment, VariantData } from '@/lib/types'
import Loading from './loading'

function AnalysisContent() {
  const searchParams = useSearchParams()
  const { user } = useAuth()
  
  const [experiments, setExperiments] = useState<Experiment[]>([])
  const [selectedExpId, setSelectedExpId] = useState<string>(searchParams.get('experiment') || '')
  const [variants, setVariants] = useState<VariantData[]>([])
  const [selectedExperiment, setSelectedExperiment] = useState<Experiment | null>(null)
  const [isLoading, setIsLoading] = useState(true)
  const [selectedVariantIndex, setSelectedVariantIndex] = useState<number | null>(null)

  // Load experiments list
  useEffect(() => {
    async function loadExperiments() {
      if (!user) return
      
      try {
        const res = await fetch(`${process.env.NEXT_PUBLIC_BACKEND_URL || 'http://localhost:8000'}/api/experiments?userId=${user.id}`, {
          credentials: 'include'
        })
        const data = await res.json()
        if (data.success) {
          const validExps = data.experiments.filter((e: Experiment) => e.validationStatus === 'valid')
          setExperiments(validExps)
          
          // Auto-select if only one or from URL param
          if (validExps.length === 1) {
            setSelectedExpId(validExps[0].id)
          } else if (searchParams.get('experiment')) {
            setSelectedExpId(searchParams.get('experiment') || '')
          }
        }
      } catch (err) {
        console.error('Failed to load experiments:', err)
      } finally {
        setIsLoading(false)
      }
    }

    loadExperiments()
  }, [user, searchParams])

  // Load variants when experiment selected
  useEffect(() => {
    async function loadVariants() {
      if (!selectedExpId) {
        setVariants([])
        setSelectedExperiment(null)
        return
      }

      try {
        const res = await fetch(`${process.env.NEXT_PUBLIC_BACKEND_URL || 'http://localhost:8000'}/api/experiments/${selectedExpId}`, {
          credentials: 'include'
        })
        const data = await res.json()
        
        if (data.success) {
          setSelectedExperiment(data.experiment)
          setVariants(data.variants || [])
        }
      } catch (err) {
        console.error('Failed to load variants:', err)
      }
    }

    loadVariants()
  }, [selectedExpId])

  if (isLoading) {
    return null
  }

  const passedVariants = variants.filter(v => v.qcStatus === 'passed')
  const topPerformers = [...passedVariants]
    .sort((a, b) => b.activityScore - a.activityScore)
    .slice(0, 10)

  return (
    <div className="space-y-6">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <h1 className="text-2xl font-bold text-foreground">Analysis Dashboard</h1>
          <p className="text-muted-foreground mt-1">
            Visualize and analyze your directed evolution experiments
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
              No experiments with uploaded data found. Create an experiment and upload data to see analysis.
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
      ) : variants.length === 0 ? (
        <Card>
          <CardContent className="py-12 text-center">
            <Dna className="h-12 w-12 mx-auto text-muted-foreground/50 mb-3" />
            <p className="text-muted-foreground">
              No variant data uploaded for this experiment. Upload data to see analysis.
            </p>
          </CardContent>
        </Card>
      ) : (
        <>
          {/* Summary stats */}
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
                <p className="text-2xl font-bold text-accent">{passedVariants.length}</p>
              </CardContent>
            </Card>
            <Card>
              <CardContent className="pt-6">
                <p className="text-sm text-muted-foreground">Avg Activity</p>
                <p className="text-2xl font-bold">
                  {(passedVariants.reduce((s, v) => s + v.activityScore, 0) / passedVariants.length || 0).toFixed(2)}
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
                  {Math.max(...passedVariants.map(v => v.activityScore), 0).toFixed(2)}
                </p>
              </CardContent>
            </Card>
          </div>

          {/* Visualizations */}
          <Tabs defaultValue="overview">
            <TabsList>
              <TabsTrigger value="overview">Overview</TabsTrigger>
              <TabsTrigger value="performers">Top Performers</TabsTrigger>
              <TabsTrigger value="mutations">Mutation Fingerprint</TabsTrigger>
              <TabsTrigger value="landscape">Activity Landscape</TabsTrigger>
            </TabsList>

            <TabsContent value="overview" className="space-y-6 mt-6">
              <ActivityDistributionChart variants={passedVariants} />
              
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
                    onSelectVariant={(idx) => setSelectedVariantIndex(idx)}
                    selectedVariant={selectedVariantIndex}
                  />
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="performers" className="mt-6">
              <Card>
                <CardHeader>
                  <CardTitle>Top Performing Variants</CardTitle>
                  <CardDescription>
                    Detailed view of the top 10 variants by unified activity score
                  </CardDescription>
                </CardHeader>
                <CardContent>
                  <TopPerformersTable 
                    variants={topPerformers} 
                    onSelectVariant={(idx) => setSelectedVariantIndex(idx)}
                    selectedVariant={selectedVariantIndex}
                    detailed
                  />
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="mutations" className="mt-6">
              <MutationFingerprint 
                variants={topPerformers}
                wtSequence={selectedExperiment?.protein?.sequence || ''}
                selectedVariantIndex={selectedVariantIndex}
                onSelectVariant={(idx) => setSelectedVariantIndex(idx)}
              />
            </TabsContent>

            <TabsContent value="landscape" className="mt-6">
              <ActivityLandscape variants={passedVariants} />
            </TabsContent>
          </Tabs>
        </>
      )}
    </div>
  )
}

export default function AnalysisPage() {
  return (
    <Suspense fallback={<Loading />}>
      <AnalysisContent />
    </Suspense>
  )
}
