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

import java.util.Random
import java.util.concurrent.{Callable, ExecutorService}
import java.util.concurrent.atomic.AtomicIntegerArray
import breeze.stats.distributions.{RandBasis, Binomial}
import org.apache.commons.math3.random.MersenneTwister

import collection.JavaConversions._

/**
 * The main class implementing C4 and ClusterWild! (both BSP and asynchronous versions).
 *
 * C4 is a serializable parallel correlation clustering algorithm.
 * Given a graph and an ordering, it ensures the same clustering is obtained as the serial KwikCluster algorithm.
 * hence, it emphasizes consistency, while also trying to maximize concurrency.
 *
 * ClusterWild! emphasizes concurrency over consistency, by performing the KwikCluster greedy peeling process in parallel without checking for potential conflicts.
 * As such, it may make mistakes and return a different outcome as KwikCluster.
 * However, we show in our paper that the expectedd error of ClusterWild! can be bounded and is close to that of serial KwikCluster.
 *
 * @param graph Graph on which to run correlation clustering.
 */
class ParallelCorrelationClustering (graph: SparseGraph) {
  private val nElements: Int = graph.numNodes()

  /**
   * Actual function that performs the parallel correlation clustering.
   * @param nThreads Number of threads to use for parallel correlation clustering.
   * @param threads A thread pool with at least nThreads threads.
   * @param ordering The order in which to process the vertices: order(i) is the ith vertex to be processed
   * @param invOrder The inverse function of order: invOrder(v) is equal to order(i) iff order(i) = v
   * @param clusterID Output array: |clusterID(v)| is the cluster to which vertex v belongs; the id is positive iff v is a cluster center.
   * @param clusterID_atomic An atomic array corresponding to clusterID (see appendix B.1 of paper).
   * @param doClusterWild Boolean indicating whether to do ClusterWild! or C4.
   * @param doBSP Boolean indicating whether to do the BSP or asynchronous version.
   * @param epsilon Determines the sampling rate per round: epsilon / maxDegree of the remaining vertices are sampled.
   * @param delta Delta is used to compute the number of rounds required for max degree to halve, with probability 1-delta.
   * @param randSeed A random seed for the random sampling process.
   * @return Tuple with number of edges, number of vertices that have to wait in C4, number of vertices that were proposed but rejected as cluster centers in C4, and number of BSP rounds.
   */
  def run(nThreads: Int,
          threads: ExecutorService,
          ordering: Array[Int],
          invOrder: Array[Int],
          clusterID: Array[Int],
          clusterID_atomic: AtomicIntegerArray,
          doClusterWild: Boolean = false,
          doBSP: Boolean = true,
          epsilon : Double = 0.5,
          delta : Double = 0.1,
          randSeed : Long = 94720)
  : (Long, Int, Int, Int) = {

    val eps = epsilon

    var maxDeg : Int = if (doBSP) maxDegreeInit(nThreads, threads) else 0
    var batchBegPoint : Int = 0 // batch begins at batchBegPoint
    var batchEndPoint : Int = 0 // batch ends   at batchEndPoint-1, i.e. *before* batchEndPoint

    var maxDegHalfLife : Int = math.ceil(2.0 / eps * math.log(nElements.toDouble * math.log(maxDeg.toDouble) / delta)).toInt
    var maxDegHalfLifeCounter : Int = 0

    val numEdges : Array[Long] = Array.fill(nThreads)(0)
    val numWaits : Array[Int]  = Array.fill(nThreads)(0)
    val numRjcts : Array[Int]  = Array.fill(nThreads)(0)

    val rand : Random = new Random(randSeed)

    var iteration = 0
    while (batchEndPoint < nElements){
      batchBegPoint = batchEndPoint
      batchEndPoint = if (maxDeg==0 || !doBSP) nElements else
        math.min(nElements, batchEndPoint + (new Binomial(nElements-batchEndPoint, eps / maxDeg.toDouble)(new RandBasis(new MersenneTwister(rand.nextInt())))).draw())
      // batchEndPoint = if (maxDeg==0 || !doBSP) nElements else
      //   math.min(nElements, batchEndPoint + (eps / maxDeg.toDouble * (nElements-batchEndPoint).toDouble).toInt)

      val tasks_Peel = (0 until nThreads).map(threadID => new Callable[Unit]{
        override def call() = {
          var index = threadID + batchBegPoint
          var v = 0
          var isNewCenter_needsWait = (false, false)
          var numEdgesLocal : Long = 0
          var numWaitsLocal : Int  = 0
          var numRjctsLocal : Int  = 0

          while (index < batchEndPoint){
            v = ordering(index)
            // We make this vertex a center if:
            // (CW) it has not been assigned to a center earlier than batchBegPoint
            // (C4) serializability check holds
            if (clusterID(v) == 0 || (doBSP && doClusterWild && clusterID(v) <= -(batchBegPoint+1))) {
              numEdgesLocal += (if (!doClusterWild) graph.nSuccessor(v) else 0)
              isNewCenter_needsWait = processVertex(v, invOrder, clusterID, clusterID_atomic, batchBegPoint, doClusterWild, doBSP)
              if (isNewCenter_needsWait._1) {
                peelCenter(v, invOrder, clusterID, clusterID_atomic)
                numEdgesLocal += graph.nSuccessor(v)
              }else{
                numRjctsLocal += 1
              }
              if (isNewCenter_needsWait._2) numWaitsLocal += 1
            }
            index += nThreads
          }
          numEdges(threadID) += numEdgesLocal
          numWaits(threadID) += numWaitsLocal
          numRjcts(threadID) += numRjctsLocal
        }
      })
      threads.invokeAll(tasks_Peel)

      maxDegHalfLifeCounter += 1
      if (doBSP && maxDegHalfLifeCounter >= maxDegHalfLife){
        maxDeg = math.ceil(maxDeg.toFloat / 2.0).toInt
        maxDegHalfLifeCounter = 0
        if (maxDeg<=1){
          maxDeg = maxDegree(nThreads, threads, clusterID, ordering, batchEndPoint)
          maxDegHalfLife = math.ceil(2.0 / eps * math.log((nElements-batchEndPoint).toDouble * math.log(maxDeg.toDouble) / delta)).toInt
        }
      }
      iteration += 1
    }

    ( numEdges.reduce(_+_), numWaits.reduce(_+_), numRjcts.reduce(_+_), iteration )
  }

  /**
   * Method encoding a transaction on each vertex v.
   * Within this function, we determine if the vertex should be made into a cluster center;.
   * @param v Vertex to process
   * @param invOrder The inverse function of order: invOrder(v) is equal to order(i) iff order(i) = v
   * @param clusterID Output cluster id of vertices.
   * @param clusterID_atomic Atomic array corresponding to clusterID (see appendix B.1).
   * @param batchBegPoint Start point of the vertex
   * @param doClusterWild Boolean indicating whether to do ClusterWild! or C4.
   * @param doBSP Boolean indicating whether to do the BSP or asynchronous version.
   * @return Tuple of booleans indicating whether the vertex is a center, and whether we had to wait (for C4) to process the vertex.
   */
  private def processVertex(v: Int,
                    invOrder: Array[Int],
                    clusterID: Array[Int],
                    clusterID_atomic: AtomicIntegerArray,
                    batchBegPoint: Int,
                    doClusterWild: Boolean,
                    doBSP: Boolean)
  : (Boolean, Boolean) = {
    // Check if *really*
    // (CW) clusterID == 0 or clusterID <= -(batchBegPoint+1)
    // (C4) clusterID == 0 and no earlier neighbor is center
    var vClusterID = clusterID_atomic.get(v)
    clusterID(v) = vClusterID
    if (vClusterID == 0 || (doBSP && doClusterWild && vClusterID <= -(batchBegPoint+1))){
      val nSucc = graph.nSuccessor(v)
      val vOrder = invOrder(v)

      var isCovered = !doClusterWild
      var needsWait = false

      if (!doClusterWild){
        var u = 0
        var ui = 0
        var uOrder = 0
        var uClusterID = 0

        ui = 0
        isCovered = false
        while (ui < nSucc && !isCovered){
          u = graph.succ(v, ui)
          uOrder = invOrder(u)
          if (uOrder < vOrder){
            // Neighbor with earlier ordering, check if it is cluster center
            uClusterID = clusterID(u)
            if (uClusterID == 0) needsWait = true
            while (uClusterID == 0){
              // Check if *really* ==0
              uClusterID = clusterID_atomic.get(u)
            }
            vClusterID = clusterID(v)
            if (vClusterID < 0){
              // Someone else has now covered this vertex
              isCovered = true
            }else if (uClusterID > 0){
              // Earlier neighbor is center!
              isCovered = true
            }
          }
          ui += 1
        }

      }

      if (!isCovered){
        clusterID(v) = vOrder+1
        clusterID_atomic.set(v, vOrder+1)
      }
      (!isCovered, needsWait)

    }else{
      (false, false)
    }
  }

  /**
   * Method for creating a cluster center at vertex v.
   * @param v New cluster center
   * @param invOrder The inverse function of order: invOrder(v) is equal to order(i) iff order(i) = v
   * @param clusterID Output cluster id of vertices.
   * @param clusterID_atomic Atomic array corresponding to clusterID (see appendix B.1).
   */
  private def peelCenter(v: Int,
                 invOrder: Array[Int],
                 clusterID: Array[Int],
                 clusterID_atomic: AtomicIntegerArray
                  ): Unit = {
    val nSucc = graph.nSuccessor(v)
    val vOrder = invOrder(v)
    val vClusterID = vOrder + 1

    var u = 0
    var ui = 0
    var uOrder = 0
    var uClusterID = 0

    while (ui < nSucc) {
      u = graph.succ(v, ui)
      uOrder = invOrder(u)
      // We assign vertex u to center v if:
      // (1 ) uOrder > vOrder AND
      // (2a) u is not a center OR
      // (2b) u has been assigned to a later center
      // Note: For CW this can result in a bad assignment if u is in the same batch as v
      //       However, this should be fixed by the thread processing u (in run and processVertex methods)
      if (uOrder > vOrder) {

        uClusterID = clusterID(u)
        if (uClusterID == 0 || vClusterID < -uClusterID) {
          uClusterID = clusterID_atomic.get(u)
          while (uClusterID == 0 || vClusterID < -uClusterID) {
            clusterID_atomic.compareAndSet(u, uClusterID, -vClusterID)
            uClusterID = clusterID_atomic.get(u)
          }
          clusterID(u) = uClusterID
        }
      }
      ui += 1
    }

  }

  /**
   * Initial computation of maximum degree of the graph.
   * @param nThreads Number of threads to use for computing maximum degree of graph.
   * @param threads Thread pool containing at least nThreads threads.
   * @return maximum degree of graph.
   */
  private def maxDegreeInit(nThreads: Int, threads: ExecutorService) : Int = {
    val tasks = (0 until nThreads).map(threadID => new Callable[Int] {
      override def call() = {
        var localMaxDeg : Int = 0
        var i = threadID
        while (i < nElements){
          localMaxDeg = if (localMaxDeg > graph.nSuccessor(i)) localMaxDeg else graph.nSuccessor(i)
          i += nThreads
        }
        localMaxDeg
      }
    })
    val maxDeg = threads.invokeAll(tasks).map(_.get()).max
    maxDeg
  }

  /**
   * Computes maximum degree of a partially processed / clustered graph.
   * @param nThreads Number of threads to use for computing maximum degree of graph.
   * @param threads Thread pool containing at least nThreads threads.
   * @param clusterID Cluster id of (partially clustered) vertices.
   * @param ordering The order in which to process the vertices: order(i) is the ith vertex to be processed.
   * @param firstUnprocessed The index of the first vertex that (may be) unclustered.
   * @return maximum degree of the unprocessed / unclustered graph.
   */
  private def maxDegree(nThreads: Int, threads: ExecutorService, clusterID: Array[Int], ordering: Array[Int], firstUnprocessed: Int) : Int = {
    val tasks = (0 until nThreads).map(threadID => new Callable[Int] {
      override def call() = {
        var localMaxDeg : Int = 0
        var index = threadID + firstUnprocessed
        var j = 0
        var nSucc = 0
        var v = 0
        var vDeg: Int = 0
        while (index < nElements){
          v = ordering(index)
          vDeg = 0
          if (clusterID(v) == 0) {
            nSucc = graph.nSuccessor(v)
            j = 0
            while (j < nSucc) {
              if (clusterID(graph.succ(v, j)) == 0) vDeg += 1
              j += 1
            }
          }
          localMaxDeg = if (localMaxDeg > vDeg) localMaxDeg else vDeg
          index += nThreads
        }
        localMaxDeg
      }
    })
    val maxDeg = threads.invokeAll(tasks).map(_.get()).max
    maxDeg
  }
}
