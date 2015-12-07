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

import scala.collection.JavaConversions._
import scala.collection.mutable.HashSet

/**
 * Parallel correlation clustering algorithm of [6] Chierichetti, Dalvi, and Kumar (2014).
 * At each BSP round, the algorithm samples a small batch (inversely proportional to the max degree of the unprocessed graph) of unprocessed vertices.
 * Sampled vertices that are not neighbors with any other sampled vertex are made into cluster centers.
 * Neighbors of the new cluster centers that are yet unassigned (not already clustered, nor made into a new cluster center) are then assigned to these new cluster centers.
 * @param graph Graph on which to run correlation clustering.
 */
class ParallelCorrelationClustering_CDK (graph: SparseGraph) {
  /**
   * Number of elements in the graph.
   */
  private val nElements: Int = graph.numNodes()

  /**
   * Actual function that performs the parallel correlation clustering.
   * @param nThreads Number of threads to use for parallel correlation clustering.
   * @param threads A thread pool with at least nThreads threads.
   * @param clusterID Output array with cluster IDs (see description of SequentialRandomGreedyPeeling).
   * @param clusterID_atomic An atomic array corresponding to clusterID (see appendix B.1 of paper).
   * @param randSeed A random seed for the random sampling process.
   * @param alwaysComputeMaxDeg Boolean indicating whether to recompute the maximum degree after each round (see appendix B.2).
   * @param epsilon Determines the sampling rate per round: epsilon / maxDegree of the remaining vertices are sampled.
   * @param delta Delta is used to compute the number of rounds required for max degree to halve, with probability 1-delta.
   * @return Number of edges examined, and number of BSP rounds.
   */
  def run(nThreads: Int,
          threads: ExecutorService,
          clusterID: Array[Int],
          clusterID_atomic: AtomicIntegerArray,
          randSeed: Long = 94720,
          alwaysComputeMaxDeg: Boolean = false,
          epsilon: Double = 0.5,
          delta: Double = 0.1)
  : (Long, Int) = {

    var eps = epsilon
    if (epsilon >= 1.0){
      eps = 0.9999
      System.err.println(s"epsilon set too high at $epsilon, resetting to epsilon = $eps")
    }

    var maxDeg: Int = maxDegreeInit(nThreads, threads)
    var maxDegHalfLife : Int = math.ceil(2.0 / eps * math.log(nElements.toDouble * math.log(maxDeg.toDouble) / delta)).toInt
    var maxDegHalfLifeCounter : Int = 0

    val numEdges : Array[Long] = Array.fill(nThreads)(0)

    var hasUnprocessed : Boolean = true

    // TODO: for performance we use Array, but lost updates may be a *minor* issue
    val isActive : Array[Boolean] = Array.fill(nElements)(false)

    val masterRand : Random = new Random(randSeed)
    val rands : Array[Random] = (0 until nThreads).map(_ => new Random(masterRand.nextLong())).toArray

    var iteration = 0

    // Create HashSets of unprocessed vertices
    val tasks_createUnprocessedHashSets = (0 until nThreads).map(threadID => new Callable[HashSet[Int]]{
      override def call() = {
        val localUnprocessed : HashSet[Int] = new HashSet[Int]()
        var index = threadID
        while (index < nElements){
          localUnprocessed += index
          index += nThreads
        }
        localUnprocessed
      }
    })
    val unprocessedVertices : Array[HashSet[Int]] = threads.invokeAll(tasks_createUnprocessedHashSets).map(_.get).toArray

    iteration = 0
    while (hasUnprocessed){
      // Create active set
      val tasks_activeSet = (0 until nThreads).map(threadID => new Callable[HashSet[Int]]{
        override def call() = {
          val localActiveSet : HashSet[Int] = unprocessedVertices(threadID).filter(
            v => (rands(threadID).nextDouble() < eps / maxDeg.toDouble && clusterID_atomic.get(v)==0))
          localActiveSet.foreach(v => isActive(v) = true)

          localActiveSet
        }
      })
      val activeSets : Array[HashSet[Int]] = threads.invokeAll(tasks_activeSet).map(_.get).toArray

      // Check for friends
      val tasks_checkFriends = (0 until nThreads).map(threadID => new Callable[HashSet[Int]] {
        override def call() = {
          activeSets(threadID).filter(v => {
            var ui = 0
            var u = 0
            val nSucc = graph.nSuccessor(v)
            var hasFriend = false
            while (ui < nSucc && !hasFriend) {
              u = graph.succ(v, ui)
              if (isActive(u)) hasFriend = true
              ui += 1
            }
            numEdges(threadID) += ui
            hasFriend
          })
        }
      })
      val friendlyVertices : Array[HashSet[Int]] = threads.invokeAll(tasks_checkFriends).map(_.get).toArray

      // Remove friendly vertices from isActive, activeSets
      val tasks_removeFriendlies = (0 until nThreads).map(threadID => new Callable[Unit] {
        override def call() = {
          activeSets(threadID) --= friendlyVertices(threadID)
          friendlyVertices(threadID).foreach(v => isActive(v) = false)
        }
      })
      threads.invokeAll(tasks_removeFriendlies)

      // Peel away active set; also clear isActive
      val tasks_peelCenters = (0 until nThreads).map(threadID => new Callable[Unit] {
        override def call() = {
          activeSets(threadID).foreach(v => {
            clusterID(v) = v+1
            clusterID_atomic.set(v, clusterID(v))
            peelCenter(v, clusterID, clusterID_atomic)
            numEdges(threadID) += graph.nSuccessor(v)
            isActive(v) = false
          })
        }
      })
      threads.invokeAll(tasks_peelCenters)

      val tasks_removeProcessed = (0 until nThreads).map(threadID => new Callable[Unit] {
        override def call() = {
          unprocessedVertices(threadID) --= unprocessedVertices(threadID).filter(v => clusterID(v) != 0)
        }
      })
      threads.invokeAll(tasks_removeProcessed)

      hasUnprocessed = unprocessedVertices.map(_.size > 0).reduce(_|_)

      maxDegHalfLifeCounter += 1
      if (maxDegHalfLifeCounter >= maxDegHalfLife){
        maxDeg = math.ceil(maxDeg.toFloat / 2.0).toInt
        maxDegHalfLifeCounter = 0
        if (maxDeg<=1){
          maxDeg = maxDegree(nThreads, threads, clusterID, unprocessedVertices)
          maxDegHalfLife = math.ceil(2.0 / eps * math.log(nElements.toDouble * math.log(maxDeg.toDouble) / delta)).toInt
        }
      }
      iteration += 1
    }

    (numEdges.sum, iteration)
  }

  /**
   * Method for creating a cluster center at vertex v.
   * @param v New cluster center
   * @param clusterID Output cluster id of vertices.
   * @param clusterID_atomic Atomic array corresponding to clusterID (see appendix B.1).
   */
  private def peelCenter(v: Int,
                 clusterID: Array[Int],
                 clusterID_atomic: AtomicIntegerArray
                  ): Unit = {
    val nSucc = graph.nSuccessor(v)
    val vOrder = v
    val vClusterID = vOrder + 1

    var u = 0
    var ui = 0
    var uOrder = 0
    var uClusterID = 0

    while (ui < nSucc) {
      u = graph.succ(v, ui)
      uOrder = u
      // We assign vertex u to center v if:
      // (1 ) uOrder > vOrder AND
      // (2a) u is not a center OR
      // (2b) u has been assigned to a later center
      // Note: CDK ensures that if v is a center, then it has no neighbors in this batch.
      //       Thus, condition (1) implies u has not been activated,
      //       and hence u must be a spoke because it is neighbor to v
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
   * @param unprocessedVertices Set of unprocessed vertices: unprocessed(threadID) is the set of vertices that thread threadID has yet to cluster.
   * @return maximum degree of the unprocessed / unclustered graph.
   */
  private def maxDegree(nThreads: Int, threads: ExecutorService, clusterID: Array[Int], unprocessedVertices: Array[HashSet[Int]]) : Int = {
    val tasks = (0 until nThreads).map(threadID => new Callable[Int] {
      override def call() = {
        var localMaxDeg : Int = 0
        var j = 0
        var nSucc = 0
        var vDeg: Int = 0
        unprocessedVertices(threadID).foreach(v => {
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
        })
        localMaxDeg
      }
    })
    val maxDeg = threads.invokeAll(tasks).map(_.get()).max
    maxDeg
  }
}
