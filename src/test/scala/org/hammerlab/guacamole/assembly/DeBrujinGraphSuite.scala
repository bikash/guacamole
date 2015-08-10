package org.hammerlab.guacamole.assembly

import org.hammerlab.guacamole.util.{ GuacFunSuite, TestUtil }
import org.scalatest.Matchers

class DeBrujinGraphSuite extends GuacFunSuite with Matchers {

  lazy val reads = TestUtil.loadReads(sc, "NA12878_S1-chr1-10000.sam").mappedReads.collect()
  lazy val smallWindowReads = reads.filter(read => read.start > 10000 && read.start < 10100)

  lazy val snpReads = TestUtil.loadReads(sc, "assemble-reads-chr2-130508072-130508320.sam").mappedReads.collect()

  lazy val smallWindowSequences = {
    smallWindowReads.map(_.sequence)
  }

  test("Kmer.mergeKmers") {
    val kmers = Seq("TTTC", "TTCC", "TCCC", "CCCC")
    val longerKmer = DeBrujinGraph.mergeKmers(kmers)

    longerKmer.length should be (7)
    longerKmer.toString should be("TTTCCCC")

  }

  test("build graph") {

    val sequence = "TCATCTCAAAAGAGATCGA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 8)

    val firstKmer = "TCATCTCA"
    val nextKmer = "CATCTCAA"
    val lastKmer = "GAGATCGA"

    assert(graph.kmerCounts.contains(firstKmer))
    assert(graph.kmerCounts.contains(nextKmer))
    assert(graph.kmerCounts.contains(lastKmer))
  }

  test("build graph test prefix/suffix") {

    val sequence = "TCATCTCAAAAGAGATCGA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 8)

    val firstKmer = "TCATCTCA"
    val nextKmer = "CATCTCAA"
    val lastKmer = "GAGATCGA"

    graph.kmerPrefix("TCATCTCA") should be("TCATCTC")
    graph.kmerPrefix("CATCTCAA") should be("CATCTCA")
    graph.kmerPrefix("GAGATCGA") should be("GAGATCG")

    graph.kmerSuffix("TCATCTCA") should be("CATCTCA")
    graph.kmerSuffix("CATCTCAA") should be("ATCTCAA")
    graph.kmerSuffix("GAGATCGA") should be("AGATCGA")

  }

  test("build graph with short kmers and correct counts") {

    val sequence = "TCATCTTAAAAGACATAAA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 3)

    val firstKmer = "TCA"
    val nextKmer = "CAT"
    val lastKmer = "AAA"

    assert(graph.kmerCounts(firstKmer) === 1)
    assert(graph.kmerCounts(nextKmer) === 2)
    assert(graph.kmerCounts(lastKmer) === 3)
  }

  test("build graph with short kmers and correct children/parents") {

    val sequence = "TCATCTTAAAAGACATAAA"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 3)

    val firstKmer = "TCA"
    val nextKmer = "CAT"
    val lastKmer = "AAA"

    val tcaChildren = graph.children("TCA")
    graph.kmerSuffix("TCA") should be("CA")
    tcaChildren.length should be(1) // CA is the suffix and CAT is the only kmer
    tcaChildren(0) should be("CAT")

    val tcaParents = graph.parents("TCA")
    graph.kmerPrefix("TCA") should be("TC")
    tcaParents.length should be(1) // TC is the prefix and ATC is the only kmer
    tcaParents(0) should be("ATC")

    val catParents = graph.parents("CAT")
    graph.kmerPrefix("CAT") should be("CA")
    catParents.length should be(2) // CA is the prefix, TCA and ACA are parents
    catParents(0) should be("TCA")
    catParents(1) should be("ACA")

    val catChildren = graph.children("CAT")
    graph.kmerSuffix("CAT") should be("AT")
    catChildren.length should be(2) // CA is the suffix, ATC and ATA are children
  }

  test("build graph with all unique kmers") {
    val sequence = "AAATCCCTTTTA"
    val kmerSize = 4
    val graph = DeBrujinGraph(Seq(sequence), kmerSize)
    graph.kmerCounts.keys.size should be(sequence.length - kmerSize + 1)

    graph.kmerCounts.foreach(_._2 should be(1))

  }

  test("find forward unique path; full graph") {

    val sequence = "AAATCCCTGGGT"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 4)

    val firstKmer = "AAAT"

    val merageableForward = graph.mergeForward(firstKmer)
    merageableForward.size should be(9)

    //TestUtil.assertBases(graph.mergeKmers(merageableForward.reverse), sequence)
    DeBrujinGraph.mergeKmers(merageableForward.reverse) should be(sequence)
  }

  test("find backward unique path; full graph") {

    val sequence = "AAATCCCTGGGT"
    val graph = DeBrujinGraph(Seq(sequence), kmerSize = 4)

    val firstKmer = "GGGT"

    val merageableReverse = graph.mergeBackward(firstKmer)
    merageableReverse.size should be(9)
    //TestUtil.assertBases(, sequence)
    DeBrujinGraph.mergeKmers(merageableReverse) should be(sequence)
  }

  test("find forward unique path; with bubble") {

    val sequence = "AAATCCCTGGGT"
    val variantSequence = "AAATCCCTGGAT"
    val graph = DeBrujinGraph(Seq(sequence, variantSequence), kmerSize = 4)

    val firstKmer = "AAAT"

    val merageableForward = graph.mergeForward(firstKmer)
    merageableForward.size should be(7)
    TestUtil.assertBases(DeBrujinGraph.mergeKmers(merageableForward.reverse), "AAATCCCTGG")
  }

  test("find forward unique path; with bubble in middle") {

    val sequence = "AAATCCCTGGGT"
    val variantSequence = "AAATCGCTGGGT"
    val graph = DeBrujinGraph(Seq(sequence, variantSequence), kmerSize = 4)

    val firstKmer = "AAAT"

    val merageableForward = graph.mergeForward(firstKmer)
    merageableForward.size should be(2)
    TestUtil.assertBases(DeBrujinGraph.mergeKmers(merageableForward.reverse), "AAATC")
  }

  test("find forward unique path; with bubble in first kmer") {

    val sequence = "AAATCCCTGGGT"
    val variantSequence = "ACATCCCTGGGT"
    val graph = DeBrujinGraph(Seq(sequence, variantSequence), kmerSize = 4)

    val firstKmer = "AAAT"

    val merageableForward = graph.mergeForward(firstKmer)
    merageableForward.size should be(2)
    TestUtil.assertBases(DeBrujinGraph.mergeKmers(merageableForward.reverse), "AAATC")
  }

  test("find backward unique path; with bubble") {

    val sequence = "AAATCCCTGGGT"
    val variantSequence = "AAATCCCTGGAT"
    val graph = DeBrujinGraph(Seq(sequence, variantSequence), kmerSize = 4)

    val seq1End = "GGGT"

    val seq1MerageableReverse = graph.mergeBackward(seq1End)
    seq1MerageableReverse.size should be(2)
    TestUtil.assertBases(DeBrujinGraph.mergeKmers(seq1MerageableReverse), "TGGGT")

    val seq2End = "GGAT"
    val seq2MerageableReverse = graph.mergeBackward(seq2End)
    seq2MerageableReverse.size should be(2)
    TestUtil.assertBases(DeBrujinGraph.mergeKmers(seq2MerageableReverse), "TGGAT")

  }

  test("find backward unique path; with bubble in middle") {

    val sequence = "AAATCCCTGGGT"
    val variantSequence = "AAATCGCTGGGT"
    val graph = DeBrujinGraph(Seq(sequence, variantSequence), kmerSize = 4)

    val firstKmer = "AAAT"

    val merageableForward = graph.mergeForward(firstKmer)
    merageableForward.size should be(2)
    TestUtil.assertBases(DeBrujinGraph.mergeKmers(merageableForward.reverse), "AAATC")
  }

  test("test merge nodes; full graph") {

    val sequence = "AAATCCCTGGGT"
    val kmerSize = 4
    val graph = DeBrujinGraph(Seq(sequence), kmerSize)

    graph.kmerCounts.keys.size should be(9)

    graph.mergeNodes()
    graph.kmerCounts.keys.size should be(1)
    TestUtil.assertBases(graph.kmerCounts.keys.head, "AAATCCCTGGGT")

  }

  test("test merge nodes; with variant") {

    val sequence = "AAATCCCTGGGT"
    val variantSequence = "AAATCCCTGGAT"
    val kmerSize = 4
    val graph = DeBrujinGraph(Seq(sequence, variantSequence), kmerSize)

    graph.kmerCounts.keys.size should be(11)

    graph.mergeNodes()
    graph.kmerCounts.keys.size should be(3)

    assert(graph.kmerCounts.contains("AAATCCCTGG"))

    assert(graph.kmerCounts.contains("TGGGT"))
    assert(graph.kmerCounts.contains("TGGAT"))

  }

  test("find single unique path in sequence") {

    val reference =
      "GAGGATCTGCCATGGCCGGGCGAGCTGGAGGAGCGAGGAGGAGGCAGGAGGA"

    val reads =
      Seq(
        reference.substring(0, 25),
        reference.substring(5, 30),
        reference.substring(7, 32),
        reference.substring(10, 35),
        reference.substring(19, 41),
        reference.substring(22, 44),
        reference.substring(25, 47),
        reference.substring(31, 52) + "TTT"
      )

    val kmerSize = 15
    val graph: DeBrujinGraph = DeBrujinGraph(
      reads,
      kmerSize,
      minOccurrence = 1,
      mergeNodes = false
    )

    val referenceKmerSource = reference.take(kmerSize)
    val referenceKmerSink = reference.takeRight(kmerSize)
    val paths = graph.depthFirstSearch(referenceKmerSource, referenceKmerSink)

    paths.length should be(1)
    DeBrujinGraph.mergeKmers(paths(0)._1.reverse) should be(reference)

    graph.mergeNodes()
    val pathsAfterMerging = graph.depthFirstSearch(referenceKmerSource, referenceKmerSink)
    pathsAfterMerging.length should be(1)
    DeBrujinGraph.mergeKmers(pathsAfterMerging(0)._1.reverse) should be(reference)

  }

  ////  sparkTest("assemble reads with snvs") {
  ////
  ////    val kmerSize = 20
  ////    val
  ////    val graph = DeBrujinGraph(
  ////      snpReads.map(_.sequence),
  ////      kmerSize,
  ////      minOccurrence = 2,
  ////      mergeNodes = true
  ////    )
  ////
  ////    assert(graph.kmerCounts(referenceRoot) > 0)
  ////    assert(graph.kmerCounts(referenceSink) > 0)
  ////
  ////    val paths = graph.depthFirstSearch(referenceRoot, referenceSink)
  ////    paths.length === 4
  ////  }
}
