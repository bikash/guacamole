package org.hammerlab.guacamole.assembly

import org.hammerlab.guacamole.reads.MDTagUtils
import org.hammerlab.guacamole.{Bases, Common, ReadSet}


object DeBrujinGraphSimpleAssembly extends App {
  
  override def main(args: Array[String]) = {
//    val kmerSize = 17
//
    println("DeBrujinGraphSimpleAssembly")
    val sc = Common.createSparkContext()

//    val reference =
//      "GAGGATCTGCCATGGCCGGGCGAGCTGGAGGAGCGAGGAGGAGGCAGGAGGA"
//
//    val sequences =
//      Seq(
//        reference.substring(0, 25),
//        reference.substring(5, 30),
//        reference.substring(7, 32),
//        reference.substring(10, 35),
//        reference.substring(19, 41),
//        reference.substring(22, 44),
//        reference.substring(25, 47),
//        reference.substring(31, 52) + "TTT"
//      )
//
//
//    val currentGraph: DeBrujinGraph = DeBrujinGraph(
//      sequences,
//      kmerSize,
//      minOccurrence = 1,
//      mergeNodes = false
//    )
//
//    val referenceKmerSource = reference.take(kmerSize)
//    assert(referenceKmerSource.size == (kmerSize))
//    val referenceKmerSink = reference.takeRight(kmerSize)
//
//    println("Kmers in Graph")
//    currentGraph.kmerCounts.keySet.foreach(println)
//    println()
//
//    println(s"Path start: $referenceKmerSource")
//    println(s"Path end: $referenceKmerSink")
//
//    val paths = currentGraph.depthFirstSearch(referenceKmerSource, referenceKmerSink)
//    println(paths.length)
//    println(s"Reconstructed children: ${DeBrujinGraph.mergeKmers(paths(0)._1.reverse)}")


    val kmerSize = 45
    val reads = ReadSet(
      sc,
      "/Users/arahuja/data/dream/synthetic.challenge.set3.normal.withMDTags.chr2.73613071.sam",
      requireMDTagsOnMappedReads = true
    ).mappedReads.sortBy(_.start).collect

    println(reads.length)
    val reference = MDTagUtils.getReference(
      reads,
      referenceStart = 73613005,
      referenceEnd = 73613151
    )

    println(s"Reference sequence: ${Bases.basesToString(reference)}")

    val currentGraph: DeBrujinGraph = DeBrujinGraph(
      reads.map(read => Bases.basesToString(read.sequence)),
      kmerSize,
      minOccurrence = 1,
      mergeNodes = true
    )

    val referenceKmerSource = Bases.basesToString(reference.take(kmerSize))
    val referenceKmerSink = Bases.basesToString(reference.takeRight(kmerSize))

    println(s"Reference start: $referenceKmerSource")
    println(s"Reference end: $referenceKmerSink")

    println(currentGraph.children(referenceKmerSource))
    //currentGraph.printGraph

    val paths = currentGraph.depthFirstSearch(
        referenceKmerSource,
        referenceKmerSink
    )

    println(paths.length)

    println()
    println("Done!")
  }

}
