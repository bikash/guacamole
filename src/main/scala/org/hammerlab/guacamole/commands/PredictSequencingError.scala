package org.hammerlab.guacamole.commands

import org.apache.spark.SparkContext
import org.apache.spark.mllib.classification.LogisticRegressionWithSGD
import org.apache.spark.mllib.evaluation.BinaryClassificationMetrics
import org.apache.spark.mllib.linalg.{ Vector, Vectors }
import org.apache.spark.mllib.regression.LabeledPoint
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.rdd.ADAMContext
import org.bdgenomics.formats.avro.Variant
import org.hammerlab.guacamole.Common.Arguments.Reads
import org.hammerlab.guacamole._
import org.hammerlab.guacamole.pileup.{ Match, Mismatch, PileupElement }
import org.hammerlab.guacamole.reads.Read
import org.kohsuke.args4j.{ Option => Opt }

case class LocusErrorVector(reference: Byte,
                            alternate: Byte,
                            qualityScore: Byte,
                            leftContext: Seq[Byte],
                            rightContext: Seq[Byte],
                            longContext: Seq[Byte],
                            isMismatch: Boolean) {

  val contextLength = math.max(leftContext.size, rightContext.size)
  val pileupPairFeatures = math.pow(LocusErrorVector.NUM_BASES, 2).toInt
  val numContextFeatures = math.pow(LocusErrorVector.NUM_BASES, contextLength).toInt
  val numFeatures =
    LocusErrorVector.NUM_BASES + //ref base
      LocusErrorVector.NUM_BASES + //alt base
      pileupPairFeatures +
      1 + //quality score
      2 * numContextFeatures +
      4 //individual longContext base counts

  lazy val baseCounts: Map[String, Double] =
    longContext
      .groupBy(v => v)
      .map(kv => (Bases.baseToString(kv._1), kv._2.size.toDouble))

  def toLabeledPoint(label: Double): LabeledPoint = {
    LabeledPoint(label = label, features = toSparseVector)
  }

  def toVWInstance(label: Int): String = {
    val instance = new StringBuilder()
    instance.append(label)
    instance.append(" | ") // add label to feature delimeter

    instance.append(s"|orig ${Bases.baseToString(reference)} |alt ${Bases.baseToString(alternate)} |qual $qualityScore ")
    instance.append(s"|lc ${Bases.basesToString(leftContext)} |rc ${Bases.basesToString(rightContext)} ")

    val baseCountFeatures = baseCounts.map(kv => s"${kv._1}:${kv._2}").mkString(", ")
    instance.append(s"|b $baseCountFeatures")

    instance.toString()
  }

  def toSparseVector: Vector = {

    val els = Seq[(Int, Double)](
      (LocusErrorVector.BASES_TO_DUMMY_MAP.getOrElse(reference, 4), 1.0),
      (LocusErrorVector.NUM_BASES + LocusErrorVector.BASES_TO_DUMMY_MAP.getOrElse(alternate, 4), 1.0),
      (2 * LocusErrorVector.NUM_BASES + LocusErrorVector.seqToEncoding((Seq(reference, alternate))), 1.0),
      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures, qualityScore.toDouble),
      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures + 1 + LocusErrorVector.seqToEncoding(leftContext), 1.0),
      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures + 1 + numContextFeatures + LocusErrorVector.seqToEncoding(rightContext), 1.0),

      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures + 1 + 2 * numContextFeatures + 0, baseCounts.getOrElse("A", 0.0)),
      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures + 1 + 2 * numContextFeatures + 1, baseCounts.getOrElse("T", 0.0)),
      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures + 1 + 2 * numContextFeatures + 2, baseCounts.getOrElse("C", 0.0)),
      (2 * LocusErrorVector.NUM_BASES + pileupPairFeatures + 1 + 2 * numContextFeatures + 3, baseCounts.getOrElse("G", 0.0))
    )

    Vectors.sparse(numFeatures, els)
  }
}

object LocusErrorVector {

  val NUM_BASES = 5
  val BASES_TO_DUMMY_MAP = Map((Bases.A -> 0), (Bases.T -> 1), (Bases.C -> 2), (Bases.G -> 3), (Bases.N -> 4))

  def apply(element: PileupElement, contextLength: Int): Option[LocusErrorVector] = {
    val withinContext = (element.locus > (element.read.start + contextLength)) && (element.locus < (element.read.end - contextLength))
    val leftContext = element.leftContext(contextLength)
    val rightContext = element.rightContext(contextLength)

    val longContext = element.leftContext(20) ++ element.rightContext(20)

    val minContextLength = math.min(leftContext.size, rightContext.size)
    element.alignment match {
      case Match(ref, qual) if minContextLength == contextLength => Some(LocusErrorVector(ref, ref, qual, leftContext, rightContext, longContext, element.isMismatch))
      case Mismatch(alt, qual, ref) if minContextLength == contextLength => Some(LocusErrorVector(ref, alt, qual, leftContext, rightContext, longContext, element.isMismatch))
      case _ => None
    }
  }

  def seqToEncoding(seq: Seq[Byte]): Int = {
    seq.zipWithIndex.map(kv => math.pow(NUM_BASES, kv._2).toInt * BASES_TO_DUMMY_MAP(kv._1)).sum
  }
}

object PredictSequencingError {

  protected class Arguments extends DistributedUtil.Arguments with Reads {
    @Opt(name = "--dbsnp-file", usage = "DB-SNP VCF File")
    var dbSNPVCFFile: String = ""

    @Opt(name = "--context-length", aliases = Array("-l"),
      usage = "Context length for mutational signature")
    var contextLength: Int = 3

    @Opt(name = "--fraction", aliases = Array("-f"),
      usage = "")
    var fraction: Int = 100

    @Opt(name = "--splits", aliases = Array("-s"),
      usage = "")
    var splits: Int = 2

    @Opt(name = "--train", usage = "", forbids = Array("--dump"))
    var train: Boolean = false

    @Opt(name = "--evaluate", usage = "", forbids = Array("--dump"))
    var evaluate: Boolean = false

    @Opt(name = "--split-percent", usage = "")
    var splitPercent: Int = 75

    @Opt(name = "--dump-vw", depends = Array("--vw-path"), usage = "")
    var dump: Boolean = false

    @Opt(name = "--vw-path", depends = Array("--dump-vw"), usage = "")
    var vwPath: String = ""
  }

  object Caller extends SparkCommand[Arguments] {
    override val name = "seqerror"
    override val description = ""

    def evaluateModel(lr: LogisticRegressionWithSGD,
                      data: RDD[LabeledPoint],
                      splits: Int,
                      splitPercent: Double = 0.7): Unit = {

      val splitRDDs = data.randomSplit(Array(splitPercent, 1 - splitPercent))
      for (split <- (0 until splits)) {

        val train = splitRDDs(0)
        println(s"Training samples: ${train.count}")
        val test = splitRDDs(1)
        println(s"Test samples: ${test.count}")

        val model = lr.run(train)

        // Clear the default threshold.
        model.clearThreshold()
        println(model.weights)

        // Compute raw scores on the test set.
        val scoreAndLabels = test.map { point =>
          val score = model.predict(point.features)
          (score, point.label)
        }

        // Get evaluation metrics.
        val metrics = new BinaryClassificationMetrics(scoreAndLabels)
        val auROC = metrics.areaUnderROC()

        println(s"FOLD: ${split}, Area under ROC = ${auROC}")
      }
    }

    override def run(args: Arguments, sc: SparkContext): Unit = {

      //      val sqlContext = new SQLContext(sc)

      val filters = Read.InputFilters(mapped = true, nonDuplicate = true, passedVendorQualityChecks = true)
      val reads = Common.loadReadsFromArguments(args, sc, filters)

      val loci = Common.loci(args, reads)
      val lociPartitions = DistributedUtil.partitionLociAccordingToArgs(
        args,
        loci,
        reads.mappedReads
      )

      val (positiveAlternates, negativeAlternates) =
        buildMismatchData(
          sc,
          reads,
          args.dbSNPVCFFile,
          lociPartitions,
          args.contextLength,
          args.fraction / 100.0f)

      val splitPercent = args.splitPercent / 100.0
      if (args.dump) {
        val alternates = positiveAlternates
          .map(_.toVWInstance(1))
          .union(negativeAlternates.map(_.toVWInstance(-1)))
        val splitRDDs = alternates.randomSplit(Array(splitPercent, 1 - splitPercent))
        val train = splitRDDs(0)
        val test = splitRDDs(1)
        train.saveAsTextFile(args.vwPath + "-train")
        test.saveAsTextFile(args.vwPath + "-test")
      }

      val alternates = positiveAlternates
        .map(_.toLabeledPoint(1))
        .union(negativeAlternates.map(_.toLabeledPoint(0)))

      val lr = new LogisticRegressionWithSGD()
      lr.optimizer.setNumIterations(10)
      if (args.evaluate) {
        evaluateModel(lr, alternates, args.splits, splitPercent = splitPercent)
      }

      if (args.train) {
        //val lr = new LogisticRegressionWithLBFGS() // LogisticRegressionWithSGD(1.0, 10, 0.01, 1.0 )
        val model = lr.run(alternates)
        // Clear the default threshold.
        model.clearThreshold()
      }
    }
  }

  def buildMismatchData(sc: SparkContext,
                        reads: ReadSet,
                        dbSNPVCFFile: String,
                        lociPartitions: LociMap[Long],
                        contextLength: Int,
                        sampleFraction: Float) = {

    val dbSNPKeyedVariants: RDD[((String, Long), Variant)] = loadDbSNPPositions(sc, dbSNPVCFFile)

    val allAlternates = DistributedUtil.pileupFlatMap[((String, Long), LocusErrorVector)](
      reads.mappedReads,
      lociPartitions,
      skipEmpty = true,
      pileup => pileup
        .elements
        .flatMap(LocusErrorVector(_, contextLength))
        .map(el => ((pileup.head.read.referenceContig, pileup.locus), el)).iterator
    )

    val allAlternatesWithDBSnp = allAlternates.leftOuterJoin(dbSNPKeyedVariants)
    val positiveAlternates = allAlternatesWithDBSnp
      .filter(pe => pe._2._2.isDefined && pe._2._1.isMismatch)
      .sample(withReplacement = false, fraction = sampleFraction)

    positiveAlternates.persist()
    val numPositiveExamples = positiveAlternates.count
    println(s"Positive examples: ${numPositiveExamples}")

    val allNegativeAlternates = allAlternatesWithDBSnp
      .filter(_._2._2.isEmpty)
      .filter(_._2._1.isMismatch)

    val totalNegativeCount = allNegativeAlternates.count
    val negativeAlternates = allNegativeAlternates
      .sample(
        withReplacement = false,
        fraction = math.min(1, numPositiveExamples.toFloat / totalNegativeCount))

    negativeAlternates.persist()
    println(s"Negative examples: ${negativeAlternates.count}")
    (positiveAlternates.map(_._2._1), negativeAlternates.map(_._2._1))
  }

  def loadDbSNPPositions(sc: SparkContext, dbSNPVCFFile: String): RDD[((String, Long), Variant)] = {

    val adamContext = new ADAMContext(sc)
    val dbSNPVariants: RDD[Variant] = adamContext
      .loadVariants(dbSNPVCFFile)

    val dbSNPKeyedVariants = dbSNPVariants
      .map(v => ((v.getContig.getContigName, v.getStart.toLong), v))
    dbSNPKeyedVariants
  }
}

