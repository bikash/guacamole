<<<<<<< HEAD
///**
// * Licensed to Big Data Genomics (BDG) under one
// * or more contributor license agreements.  See the NOTICE file
// * distributed with this work for additional information
// * regarding copyright ownership.  The BDG licenses this file
// * to you under the Apache License, Version 2.0 (the
// * "License"); you may not use this file except in compliance
// * with the License.  You may obtain a copy of the License at
// *
// *     http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// */
//
//package org.hammerlab.guacamole.commands
//
//import org.apache.spark.SparkContext
//import org.apache.spark.rdd.RDD
//import org.hammerlab.guacamole.Common.Arguments.GermlineCallerArgs
//import org.hammerlab.guacamole.DistributedUtil.PerSample
//import org.hammerlab.guacamole._
//import org.hammerlab.guacamole.alignment.AlignmentState.AlignmentState
//import org.hammerlab.guacamole.alignment.{AffineGapPenaltyAlignment, AlignmentState, ReadAlignment}
//import org.hammerlab.guacamole.assembly.DeBrujinGraph
//import org.hammerlab.guacamole.filters.GenotypeFilter.GenotypeFilterArguments
//import org.hammerlab.guacamole.filters._
//import org.hammerlab.guacamole.reads.{MDTagUtils, MappedRead, Read}
//import org.hammerlab.guacamole.variants.{Allele, AlleleConversions, AlleleEvidence, CalledAllele}
//import org.hammerlab.guacamole.windowing.SlidingWindow
//import org.kohsuke.args4j.{Option => Args4jOption}
//
///**
// * *
// * Simple assembly based germline variant caller
// *
// */
//object GermlineAssemblyCaller {
//
//  protected class Arguments extends GermlineCallerArgs with GenotypeFilterArguments {
//
//    @Args4jOption(name = "--kmer-size", usage = "Length of kmer used for DeBrujin Graph assembly")
//    var kmerSize: Int = 45
//
//    @Args4jOption(name = "--snv-window-range", usage = "Number of bases before and after to check for additional matches or deletions")
//    var snvWindowRange: Int = 20
//
//    @Args4jOption(name = "--min-average-base-quality", usage = "Minimum average of base qualities in the read")
//    var minAverageBaseQuality: Int = 20
//
//    @Args4jOption(name = "--min-alignment-quality", usage = "Minimum alignment qualities of the read")
//    var minAlignmentQuality: Int = 30
//
//  }
//
//  object Caller extends SparkCommand[Arguments] {
//    override val name = "germline-assembly"
//    override val description = "call germline variants by assembling the surrounding region of reads"
//
=======
/**
 * Licensed to Big Data Genomics (BDG) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The BDG licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import org.hammerlab.guacamole.Common.Arguments.GermlineCallerArgs
import org.hammerlab.guacamole.DistributedUtil.PerSample
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.assembly.DeBrujinGraph
import org.hammerlab.guacamole.filters.GenotypeFilter.GenotypeFilterArguments
import org.hammerlab.guacamole.filters._
import org.hammerlab.guacamole.reads.{MDTagUtils, MappedRead, Read}
import org.hammerlab.guacamole.variants.{AlleleConversions, CalledAllele}
import org.hammerlab.guacamole.windowing.SlidingWindow
import org.kohsuke.args4j.{Option => Args4jOption}

/**
 * *
 * Simple assembly based germline variant caller
 *
 */
object GermlineAssemblyCaller {

  protected class Arguments extends GermlineCallerArgs with GenotypeFilterArguments {

    @Args4jOption(name = "--kmer-size", usage = "Length of kmer used for DeBrujin Graph assembly")
    var kmerSize: Int = 45

    @Args4jOption(name = "--snv-window-range", usage = "Number of bases before and after to check for additional matches or deletions")
    var snvWindowRange: Int = 20

    @Args4jOption(name = "--min-average-base-quality", usage = "Minimum average of base qualities in the read")
    var minAverageBaseQuality: Int = 20

    @Args4jOption(name = "--min-alignment-quality", usage = "Minimum alignment qualities of the read")
    var minAlignmentQuality: Int = 30

  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "germline-assembly"
    override val description = "call germline variants by assembling the surrounding region of reads"

>>>>>>> 7efc3e6... filling in germline assembly caller
//    def getVariantAlignments(alignment: ReadAlignment): Seq[(AlignmentState, Int)] = {
//      alignment
//        .alignments
//        .zipWithIndex
//        .filter(tup => AlignmentState.isGapAlignment(tup._1) || tup._1 == AlignmentState.Mismatch)
//    }
<<<<<<< HEAD
//
//    def discoverHaplotypes(graph: Option[DeBrujinGraph],
//                           windows: PerSample[SlidingWindow[MappedRead]],
//                           kmerSize: Int,
//                           minAlignmentQuality: Int,
//                           minAverageBaseQuality: Int,
//                           minOccurrence: Int = 2,
//                           expectedPloidy: Int = 2): (Option[DeBrujinGraph], Iterator[CalledAllele]) = {
//
//      val currentWindow = windows(0)
//      val locus = currentWindow.currentLocus
//
//      val newReads = currentWindow.currentRegions()
//      val filteredReads = newReads.filter(_.alignmentQuality > minAlignmentQuality)
//
//      // TODO: Should update graph instead of rebuilding it
//      // Need to keep track of reads removed from last update and reads added
//      val currentGraph: DeBrujinGraph = DeBrujinGraph(
//        filteredReads.map(read => Bases.st_.sequence),
//        kmerSize,
//        minOccurrence
//      )
//
//      val referenceStart = currentWindow.currentLocus - currentWindow.halfWindowSize
//      val referenceEnd = currentWindow.currentLocus + currentWindow.halfWindowSize
//
//      val currentReference: Array[Byte] =
//        MDTagUtils.getReference(
//          currentWindow.currentRegions(),
//          referenceStart,
//          referenceEnd
//        )
//
=======

    def discoverHaplotypes(graph: Option[DeBrujinGraph],
                           windows: PerSample[SlidingWindow[MappedRead]],
                           kmerSize: Int,
                           minAlignmentQuality: Int,
                           minAverageBaseQuality: Int,
                           minOccurrence: Int = 2,
                           expectedPloidy: Int = 2): (Option[DeBrujinGraph], Iterator[CalledAllele]) = {

      val currentWindow = windows(0)
      val locus = currentWindow.currentLocus

      val newReads = currentWindow.currentRegions()
      val filteredReads = newReads.filter(_.alignmentQuality > minAlignmentQuality)

      // Should update graph instead of rebuilding it
      // Need to keep track of reads removed from last update and reads added
      val currentGraph: DeBrujinGraph = DeBrujinGraph(
        filteredReads.map(_.sequence), 
        kmerSize, 
        minOccurrence
      )
      
      val referenceStart = currentWindow.currentLocus - currentWindow.halfWindowSize
      val referenceEnd = currentWindow.currentLocus + currentWindow.halfWindowSize
      
      val currentReference: Array[Byte] = MDTagUtils.getReference(currentWindow.currentRegions(), referenceStart, referenceEnd)

>>>>>>> 7efc3e6... filling in germline assembly caller
//      def buildVariant(locus: Long, tuple: (AlignmentState, Int)): CalledAllele = {
//
//        val variantType = tuple._1
//        val offset = tuple._2
//        val allele = Allele(Seq(currentReference((locus + offset).toInt)), Bases.stringToBases("<ALT>"))
//        val depth = filteredReads.length
//        CalledAllele(
//          currentWindow.newRegions.head.sampleName,
//          currentWindow.newRegions.head.referenceContig,
//          locus,
//          allele,
<<<<<<< HEAD
//          AlleleEvidence
//            (1,
//                depth,
//                depth,
//                depth,
//                depth,
//                filteredReads.map(_.alignmentQuality).sum / depth.toFloat, 30)
//        )
//      }
//
//      val referenceKmerSource = currentReference.take(kmerSize)
//      val referenceKmerSink = currentReference.takeRight(kmerSize)
//
//      val paths = currentGraph.depthFirstSearch(
//        referenceKmerSource,
//        referenceKmerSink
//      ).sortBy(_._2) // sort paths by score
//
//      // Take ploidy paths
//      val topPaths = paths.take(expectedPloidy).map(path => currentGraph.mergeKmers(path._1))
////      // align paths to reference
//      val alignments = topPaths.map(AffineGapPenaltyAlignment.align(_, currentReference))
////      // output variants
//      val variantAlignments = alignments.flatMap(getVariantAlignments(_).map(tup => buildVariant(locus, tup)))
//      (graph, variantAlignments.iterator)
//
//      (None, Iterator.empty)
//    }
//
//    override def run(args: Arguments, sc: SparkContext): Unit = {
//
//      val readSet = Common.loadReadsFromArguments(
//        args,
//        sc,
//        Read.InputFilters(mapped = true, nonDuplicate = true))
//
//      Common.progress(
//        "Loaded %,d mapped non-duplicate reads into %,d partitions.".format(readSet.mappedReads.count, readSet.mappedReads.partitions.length))
//
//      val loci = Common.loci(args, readSet)
//      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(args, loci, readSet.mappedReads)
//
//      val kmerSize = args.kmerSize
//      val minAlignmentQuality = args.minAlignmentQuality
//      val minAverageBaseQuality = args.minAverageBaseQuality
//
//      val lociOfInterest = DistributedUtil.pileupFlatMap[VariantLocus](
//        readSet.mappedReads,
//        lociPartitions,
//        skipEmpty = true,
//        pileup => VariantLocus(pileup).iterator
//      )
//
//      val lociOfInterestSet = LociSet.union(
//        lociOfInterest.map(
//          locus => LociSet(locus.contig, locus.locus, locus.locus + 1))
//          .collect(): _*)
//
//      val lociOfInterestPartitions = DistributedUtil.partitionLociUniformly(
//        args.parallelism,
//        lociOfInterestSet
//      )
//
//      val genotypes: RDD[CalledAllele] = DistributedUtil.windowFlatMapWithState[MappedRead, CalledAllele, Option[DeBrujinGraph]](
//        Seq(readSet.mappedReads),
//        lociPartitions,
//        skipEmpty = true,
//        halfWindowSize = args.snvWindowRange,
//        initialState = None,
//        (graph, window) => {
//          discoverHaplotypes(
//            graph,
//            window,
//            kmerSize,
//            minAlignmentQuality,
//            minAverageBaseQuality)
//        }
//      )
//
//      readSet.mappedReads.unpersist()
//
//      val filteredGenotypes = GenotypeFilter(genotypes, args).flatMap(AlleleConversions.calledAlleleToADAMGenotype)
//      Common.writeVariantsFromArguments(args, filteredGenotypes)
//      DelayedMessages.default.print()
//    }
//
//  }
//}
=======
//          AlleleEvidence(1, depth, depth, depth, depth, filteredReads.map(_.alignmentQuality).sum / depth.toFloat, 30)
//        )
//      }

      val referenceKmerSource = currentGraph.kmerPrefix(currentReference)
      val referenceKmerSink = currentGraph.kmerSuffix(currentReference)

      val paths = currentGraph.depthFirstSearch(
        referenceKmerSource,
        referenceKmerSink
      ).sortBy(_._2) // sort paths by score

      // Take ploidy paths
      val topPaths = paths.take(expectedPloidy).map(path => currentGraph.mergeKmers(path._1))
//      // align paths to reference
//      val alignments = topPaths.map(HMMAligner.align(_, currentReference))
//      // output variants
//      val variantAlignments = alignments.flatMap(getVariantAlignments(_).map(tup => buildVariant(locus, tup)))
//      (graph, variantAlignments.iterator)
      
      (None, Iterator.empty)
    }

    override def run(args: Arguments, sc: SparkContext): Unit = {

      val readSet = Common.loadReadsFromArguments(
        args,
        sc,
        Read.InputFilters(mapped = true, nonDuplicate = true))
      
      Common.progress(
        "Loaded %,d mapped non-duplicate reads into %,d partitions.".format(readSet.mappedReads.count, readSet.mappedReads.partitions.length))

      val loci = Common.loci(args, readSet)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(args, loci, readSet.mappedReads)

      val kmerSize = args.kmerSize
      val minAlignmentQuality = args.minAlignmentQuality
      val minAverageBaseQuality = args.minAverageBaseQuality
      
      val lociOfInterest = DistributedUtil.pileupFlatMap[VariantLocus](
        readSet.mappedReads,
        lociPartitions,
        skipEmpty = true,
        pileup => VariantLocus(pileup).iterator
      )
      
      val lociOfInterestSet = LociSet.union(
        lociOfInterest.map(
          locus => LociSet(locus.contig, locus.locus, locus.locus + 1))
          .collect(): _*)
      
      val lociOfInterestPartitions = DistributedUtil.partitionLociUniformly(
        args.parallelism,
        lociOfInterestSet
      )
      
      val genotypes: RDD[CalledAllele] = DistributedUtil.windowFlatMapWithState[MappedRead, CalledAllele, Option[DeBrujinGraph]](
        Seq(readSet.mappedReads),
        lociPartitions,
        skipEmpty = true,
        halfWindowSize = args.snvWindowRange,
        initialState = None,
        (graph, window) => {
          discoverHaplotypes(
            graph,
            window,
            kmerSize,
            minAlignmentQuality,
            minAverageBaseQuality)
        }
      )

      readSet.mappedReads.unpersist()

      val filteredGenotypes = GenotypeFilter(genotypes, args).flatMap(AlleleConversions.calledAlleleToADAMGenotype)
      Common.writeVariantsFromArguments(args, filteredGenotypes)
      DelayedMessages.default.print()
    }

  }
}
>>>>>>> 7efc3e6... filling in germline assembly caller
