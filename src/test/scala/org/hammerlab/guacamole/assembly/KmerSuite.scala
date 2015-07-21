//package org.hammerlab.guacamole.assembly
//
//import org.hammerlab.guacamole.Bases
//import org.hammerlab.guacamole.util.GuacFunSuite
//import org.scalatest.Matchers
//
//class KmerSuite extends GuacFunSuite with Matchers {
//
//  test("kmer last base") {
//    val tcga = Kmer("TCGA")
//    tcga.lastBase should be(Bases.A)
//
//    val acgt = Kmer("ACGT")
//    acgt.lastBase should be(Bases.T)
//  }
//
//  test("kmer tail") {
//    val tcga = Kmer("TCGA")
//    tcga.tail should be(Kmer("CGA"))
//
//    val agct = Kmer("AGCT")
//    agct.tail should be(Kmer("GCT"))
//  }
//
//  test("prefix") {
//    val tcga = Kmer("TCGA")
//    tcga.prefix should be(Bases.stringToBases("TCG"))
//
//    val acgt = Kmer("ACGT")
//    acgt.prefix should be(Bases.stringToBases(("ACG")))
//  }
//
//  test("kmer prefix") {
//    val tcga = Kmer("TCGA")
//    tcga.prefixKmer should be(Kmer("TCG"))
//
//    val acgt = Kmer("ACGT")
//    acgt.prefixKmer should be(Kmer("ACG"))
//  }
//
//  test("Kmer.mergeKmers") {
//    val kmers = Seq(Kmer("TTTC"), Kmer("TTCC"), Kmer("TCCC"), Kmer("CCCC"))
//    val longerKmer = Kmer.mergeKmers(kmers)
//
//    // longerKmer.kmerSize should be (7)
//    longerKmer.toString should be("TTTCCCC")
//
//  }
//
//}
