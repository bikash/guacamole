//package org.hammerlab.guacamole.assembly
//
//import org.hammerlab.guacamole.Bases
//
//case class Kmer(bases: Seq[Byte]) {
//
//  val kmerSize = bases.length
//
//  def prefix: Seq[Byte] = bases.slice(0, bases.length - 1)
//  def prefixKmer = Kmer(prefix)
//
//  def tail = Kmer(bases.tail)
//
//  def tailBases = bases.tail
//
//  def lastBase = bases.last
//
//  def possibleNext: Seq[Kmer] = {
//    Seq(
//      Kmer(bases.tail :+ Bases.A),
//      Kmer(bases.tail :+ Bases.T),
//      Kmer(bases.tail :+ Bases.C),
//      Kmer(bases.tail :+ Bases.G)
//    )
//  }
//
//  def possiblePrevious: Seq[Kmer] = {
//    Seq(
//      Kmer(Seq(Bases.A) ++ prefix),
//      Kmer(Seq(Bases.T) ++ prefix),
//      Kmer(Seq(Bases.C) ++ prefix),
//      Kmer(Seq(Bases.G) ++ prefix)
//    )
//  }
//
//  def :+(base: Byte): Kmer = {
//    Kmer(bases :+ base)
//  }
//
//  override def toString: String = Bases.basesToString(bases)
//}
//
//object Kmer {
//  def apply(seq: String): Kmer = {
//    Kmer(Bases.stringToBases(seq).toArray)
//  }
//
//  def mergeKmers(kmers: Seq[Kmer]): Kmer = {
//    val head = kmers.headOption.map(_.prefix).getOrElse(Seq.empty)
//    val rest = kmers.map(_.lastBase)
//    val seq = head ++ rest
//    Kmer(seq.toSeq)
//  }
//}
