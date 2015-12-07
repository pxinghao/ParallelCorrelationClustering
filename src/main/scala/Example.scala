/*
* Copyright 2015 [See AUTHORS file for list of authors]
*
*    Licensed under the Apache License, Version 2.0 (the "License");
*    you may not use this file except in compliance with the License.
*    You may obtain a copy of the License at
*
*        http://www.apache.org/licenses/LICENSE-2.0
*
*    Unless required by applicable law or agreed to in writing, software
*    distributed under the License is distributed on an "AS IS" BASIS,
*    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*    See the License for the specific language governing permissions and
*    limitations under the License.
*/

import java.util.concurrent.atomic.AtomicIntegerArray
import java.util.concurrent.Executors

import it.unimi.dsi.webgraph.BVGraph

import scala.util.Random


/**
 * Example program that runs KwikCluster, C4 (BSP, async), ClusterWild! (BSP, async), and CDK on an input graph.
 */
object Example {
  /**
   * Main program
   * @param args Set of optional command line arguments
   *
   *             inputFile:    pointer to the input graph (with no file extensions)
   *
   *             randSeedOrd:  random seed for generating random permutation
   *
   *             randSeedCDK:  random seed for CDK
   *
   *             randSeedPCC:  random seed to pass to parallel correlation clustering algorithms C4 and ClusterWild!
   *
   *             nTrials:      number of trials to run
   *
   *             maxnNThreads: maximum number of threads to test with
   *
   *             doSer:        Boolean indicating whether to do the serial KwikCluster algorithm
   *
   *             doC4:         Boolean indicating whether to do C4
   *
   *             doCW:         Boolean indicating whether to do ClusterWild!
   *
   *             doCDK:        Boolean indicating whether to do CDK multiple times
   *
   *             doCDKOnce:    Boolean indicating whether to do CDK only once, with maxNThreads
   *
   *             doBSP:        Boolean indicating whether to do BSP for C4, ClusterWild!
   *
   *             doAsync:      Boolean indicating whether to do asynchronous for C4, ClusterWild!
   *
   *             epsilon:      Determines the sampling rate per round: epsilon / maxDegree of the remaining vertices are sampled.
   *
   *             delta:        Delta is used to compute the number of rounds required for max degree to halve, with probability 1-delta.
   */
  def main(args: Array[String]) {
    val argmap = args.map { a =>
      val argPair = a.split("=")
      val name = argPair(0).toLowerCase
      val value = argPair(1)
      (name, value)
    }.toMap

    // Default options
    val inputFile = argmap.getOrElse("inputfile", "data/eswiki-2013_rand_symm")
    val randSeedOrd = argmap.getOrElse("randseedord", "94720").toInt
    val randSeedCDK = argmap.getOrElse("randseedcdk", "47").toInt
    val randSeedPCC = argmap.getOrElse("randseedpcc", "53").toInt
    val nTrials = argmap.getOrElse("ntrials", "10").toInt
    val maxNThreads = argmap.getOrElse("maxnthreads", "8").toInt
    val doSer = argmap.getOrElse("doser", "true").toBoolean
    val doC4 = argmap.getOrElse("doc4", "true").toBoolean
    val doCW = argmap.getOrElse("docw", "true").toBoolean
    val doCDK = argmap.getOrElse("docdk", "false").toBoolean
    val doCDKOnce = argmap.getOrElse("docdkonce", "true").toBoolean
    val doBSP = argmap.getOrElse("dobsp", "true").toBoolean
    val doAsync = argmap.getOrElse("doasync", "true").toBoolean
    val epsilon = argmap.getOrElse("epsilon", "0.5").toDouble
    val delta = argmap.getOrElse("delta", "0.1").toDouble

    // Create graph
    val graph: SparseGraph = new SparseGraph()
    graph.loadFromBVGraph(BVGraph.load(inputFile))
    val nVertices = graph.numNodes()
    println(s"Graph has ${graph.numArcs()} edges, ${graph.numNodes()} vertices")
    println(s"epsilon = $epsilon")

    val orderingRand : Random = new Random(randSeedOrd)
    val CDKRand : Random = new Random(randSeedCDK)
    val PCCRand : Random = new Random(randSeedPCC)

    // Cluster ids
    val seqClusterID_int   : Array[Int]         = Array.fill(nVertices)(0)
    val seqClusterID_atomic: AtomicIntegerArray = new AtomicIntegerArray(nVertices)
    val c4aClusterID_int   : Array[Int]         = Array.fill(nVertices)(0)
    val c4aClusterID_atomic: AtomicIntegerArray = new AtomicIntegerArray(nVertices)
    val c4bClusterID_int   : Array[Int]         = Array.fill(nVertices)(0)
    val c4bClusterID_atomic: AtomicIntegerArray = new AtomicIntegerArray(nVertices)
    val cwaClusterID_int   : Array[Int]         = Array.fill(nVertices)(0)
    val cwaClusterID_atomic: AtomicIntegerArray = new AtomicIntegerArray(nVertices)
    val cwbClusterID_int   : Array[Int]         = Array.fill(nVertices)(0)
    val cwbClusterID_atomic: AtomicIntegerArray = new AtomicIntegerArray(nVertices)
    val cdkClusterID_int   : Array[Int]         = Array.fill(nVertices)(0)
    val cdkClusterID_atomic: AtomicIntegerArray = new AtomicIntegerArray(nVertices)

    // Sequential statistics
    var seqStartTime: Long = 0
    var seqEndTime: Long = 0
    var seqEdgesExamined: Long = 0
    var seqObjVal: Long = 0

    // C4BS statistics
    var c4bStartTime: Long = 0
    var c4bEndTime: Long = 0
    var c4bEdgesExamined: Long = 0
    var c4bNumWaits: Int = 0
    var c4bNumRejected: Int = 0
    var c4bNumIterations: Int = 0
    var c4bObjVal: Long = 0

    // C4AS statistics
    var c4aStartTime: Long = 0
    var c4aEndTime: Long = 0
    var c4aEdgesExamined: Long = 0
    var c4aNumWaits: Int = 0
    var c4aNumRejected: Int = 0
    var c4aObjVal: Long = 0

    // CWBS statistics
    var cwbStartTime: Long = 0
    var cwbEndTime: Long = 0
    var cwbEdgesExamined: Long = 0
    var cwbNumIterations: Int = 0
    var cwbObjVal: Long = 0

    // CWAS statistics
    var cwaStartTime: Long = 0
    var cwaEndTime: Long = 0
    var cwaEdgesExamined: Long = 0
    var cwaObjVal: Long = 0

    // CDK statistics
    var cdkStartTime: Long = 0
    var cdkEndTime: Long = 0
    var cdkEdgesExamined: Long = 0
    var cdkNumIterations: Int = 0
    var cdkObjVal: Long = 0

    (0 until nTrials).foreach ( iter => {

      // Ordering
      val orderingThreads = Executors.newFixedThreadPool(maxNThreads)
      val orderingRandSeed = orderingRand.nextInt()
      val (ordering, invOrder) = AuxiliaryFunctions.parallelRandomPermutation(maxNThreads, orderingThreads, orderingRandSeed, nVertices)

      // PCC rand
      val pccRandSeed = PCCRand.nextLong()

      val threads = Executors.newFixedThreadPool(maxNThreads)

      (1 to maxNThreads).foreach(nThreads => {

        // Run sequential
        if (nThreads == 1 && doSer) {
          (0 until 1).foreach { _ => {
            System.gc()
            (0 until nVertices).foreach(v => {
              seqClusterID_int(v) = 0
              seqClusterID_atomic.set(v, 0)
            })
            seqStartTime = System.currentTimeMillis()
            val sqRunner: SequentialRandomGreedyPeeling = new SequentialRandomGreedyPeeling(graph)
            val sqResults = sqRunner.run(1, null, ordering, invOrder, seqClusterID_int)
            seqEndTime = System.currentTimeMillis()
            seqEdgesExamined = sqResults._1
          }
          }
          seqObjVal = AuxiliaryFunctions.computeObjective(maxNThreads, threads, graph, seqClusterID_int, seqClusterID_atomic)
          println(s"[Trial $iter] Sequential KwikCluster with 1 threads:")
          println("\tTime:            " + s"${seqEndTime - seqStartTime}")
          println("\t#edges examined: " + s"$seqEdgesExamined")
          println("\tObjective:       " + s"$seqObjVal")
        }

        if (doC4 && doBSP) {
          // Run C4BS
          System.gc()
          (0 until nVertices).foreach(v => {
            c4bClusterID_atomic.set(v, 0)
            c4bClusterID_int(v) = 0
          })
          c4bStartTime = System.currentTimeMillis()
          val c4bRunner: ParallelCorrelationClustering = new ParallelCorrelationClustering(graph)
          val c4bStats = c4bRunner.run(nThreads, threads, ordering, invOrder, c4bClusterID_int, c4bClusterID_atomic, doClusterWild = false, doBSP = true, epsilon = epsilon, delta = delta, randSeed = pccRandSeed)
          c4bEndTime = System.currentTimeMillis()
          c4bEdgesExamined = c4bStats._1
          c4bNumWaits = c4bStats._2
          c4bNumRejected = c4bStats._3
          c4bNumIterations = c4bStats._4
          c4bObjVal = AuxiliaryFunctions.computeObjective(maxNThreads, threads, graph, c4bClusterID_int, c4bClusterID_atomic)
          println(s"[Trial $iter] C4 (BSP) with $nThreads threads:")
          println("\tTime:            " + s"${c4bEndTime - c4bStartTime}")
          println("\t#edges examined: " + s"$c4bEdgesExamined")
          println("\tObjective:       " + s"$c4bObjVal")
          println("\t#BSP rounds:     " + s"$c4bNumIterations")
        }

        if (doC4 && doAsync) {
          // Run C4AS
          System.gc()
          (0 until nVertices).foreach(v => {
            c4aClusterID_atomic.set(v, 0)
            c4aClusterID_int(v) = 0
          })
          c4aStartTime = System.currentTimeMillis()
          val c4aRunner: ParallelCorrelationClustering = new ParallelCorrelationClustering(graph)
          val c4aStats = c4aRunner.run(nThreads, threads, ordering, invOrder, c4aClusterID_int, c4aClusterID_atomic, doClusterWild = false, doBSP = false, epsilon = epsilon, randSeed = pccRandSeed)
          c4aEndTime = System.currentTimeMillis()
          c4aEdgesExamined = c4aStats._1
          c4aNumWaits = c4aStats._2
          c4aNumRejected = c4aStats._3
          c4aObjVal = AuxiliaryFunctions.computeObjective(maxNThreads, threads, graph, c4aClusterID_int, c4aClusterID_atomic)
          println(s"[Trial $iter] C4 (async) with $nThreads threads:")
          println("\tTime:            " + s"${c4aEndTime - c4aStartTime}")
          println("\t#edges examined: " + s"$c4aEdgesExamined")
          println("\tObjective:       " + s"$c4aObjVal")
        }

        if (doCW && doBSP) {
          // Run CWBS
          System.gc()
          (0 until nVertices).foreach(v => {
            cwbClusterID_atomic.set(v, 0)
            cwbClusterID_int(v) = 0
          })
          cwbStartTime = System.currentTimeMillis()
          val cwbRunner: ParallelCorrelationClustering = new ParallelCorrelationClustering(graph)
          val cwbStats = cwbRunner.run(nThreads, threads, ordering, invOrder, cwbClusterID_int, cwbClusterID_atomic, doClusterWild = true, doBSP = true, epsilon = epsilon, delta = delta, randSeed = pccRandSeed)
          cwbEndTime = System.currentTimeMillis()
          cwbEdgesExamined = cwbStats._1
          cwbNumIterations = cwbStats._4
          cwbObjVal = AuxiliaryFunctions.computeObjective(maxNThreads, threads, graph, cwbClusterID_int, cwbClusterID_atomic)
          println(s"[Trial $iter] ClusterWild! (BSP) with $nThreads threads:")
          println("\tTime:            " + s"${cwbEndTime - cwbStartTime}")
          println("\t#edges examined: " + s"$cwbEdgesExamined")
          println("\tObjective:       " + s"$cwbObjVal")
          println("\t#BSP rounds:     " + s"$cwbNumIterations")
        }

        if (doCW && doAsync) {
          // Run CWAS
          System.gc()
          (0 until nVertices).foreach(v => {
            cwaClusterID_atomic.set(v, 0)
            cwaClusterID_int(v) = 0
          })
          cwaStartTime = System.currentTimeMillis()
          val cwaRunner: ParallelCorrelationClustering = new ParallelCorrelationClustering(graph)
          val cwaStats = cwaRunner.run(nThreads, threads, ordering, invOrder, cwaClusterID_int, cwaClusterID_atomic, doClusterWild = true, doBSP = false, epsilon = epsilon, delta = delta, randSeed = pccRandSeed)
          cwaEndTime = System.currentTimeMillis()
          cwaEdgesExamined = cwaStats._1
          cwaObjVal = AuxiliaryFunctions.computeObjective(maxNThreads, threads, graph, cwaClusterID_int, cwaClusterID_atomic)
          println(s"[Trial $iter] ClusterWild! (async) with $nThreads threads:")
          println("\tTime:            " + s"${cwaEndTime - cwaStartTime}")
          println("\t#edges examined: " + s"$cwaEdgesExamined")
          println("\tObjective:       " + s"$cwaObjVal")
        }

        if (doCDK || (doCDKOnce && nThreads == maxNThreads)) {
          // Run CDK
          System.gc()
          (0 until nVertices).foreach(v => {
            cdkClusterID_atomic.set(v, 0)
            cdkClusterID_int(v) = 0
          })
          cdkStartTime = System.currentTimeMillis()
          val cdkRunner: ParallelCorrelationClustering_CDK = new ParallelCorrelationClustering_CDK(graph)
          val cdkStats = cdkRunner.run(nThreads, threads, cdkClusterID_int, cdkClusterID_atomic, epsilon = epsilon, delta = delta, randSeed = CDKRand.nextLong())
          cdkEndTime = System.currentTimeMillis()
          cdkEdgesExamined = cdkStats._1
          cdkNumIterations = cdkStats._2
          cdkObjVal = AuxiliaryFunctions.computeObjective(maxNThreads, threads, graph, cdkClusterID_int, cdkClusterID_atomic)
          println(s"[Trial $iter] CDK with $nThreads threads:")
          println("\tTime:            " + s"${cdkEndTime - cdkStartTime}")
          println("\t#edges examined: " + s"$cdkEdgesExamined")
          println("\tObjective:       " + s"$cdkObjVal")
          println("\t#BSP rounds:     " + s"$cdkNumIterations")
        }

        //        println(s"" +
        //          s"${nThreads}\t" +
        //          s"${seqEndTime - seqStartTime}\t" +
        //          s"${seqEdgesExamined}\t" +
        //          s"${seqObjVal}\t" +
        //          s"${c4bEndTime - c4bStartTime}\t" +
        //          s"${c4bEdgesExamined}\t" +
        //          s"${c4bNumWaits}\t" +
        //          s"${c4bNumRejected}\t" +
        //          s"${c4bNumIterations}\t" +
        //          s"${c4aEndTime - c4aStartTime}\t" +
        //          s"${c4aEdgesExamined}\t" +
        //          s"${c4aNumWaits}\t" +
        //          s"${c4aNumRejected}\t" +
        //          s"${cwbEndTime - cwbStartTime}\t" +
        //          s"${cwbEdgesExamined}\t" +
        //          s"${cwbNumIterations}\t" +
        //          s"${cwbObjVal}\t" +
        //          s"${cwaEndTime - cwaStartTime}\t" +
        //          s"${cwaEdgesExamined}\t" +
        //          s"${cwaObjVal}\t" +
        //          s"${cdkEndTime - cdkStartTime}\t" +
        //          s"${cdkEdgesExamined}\t" +
        //          s"${cdkNumIterations}\t" +
        //          s"${cdkObjVal}\t" +
        //          s"")
      })


      threads.shutdown()
      orderingThreads.shutdown()

    })
  }
}