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

import java.lang.Math.abs
import java.util.Random
import java.util.concurrent.atomic.{AtomicInteger, AtomicIntegerArray}
import java.util.concurrent.{Callable, ExecutorService, Executors}

import scala.collection.JavaConversions._
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer

/**
 * Auxiliary helper functions.
 */
object AuxiliaryFunctions{
  /**
   * Randomly permute an array using the Knuth permutation algorithm.
   * @param a Array to permute.
   * @param rand Random object to use for the permutation.
   * @return Permutated array.
   */
  def KnuthPermute(a: Array[Int], rand: Random) : Array[Int] = {
    val numElements = a.length
    var i = 0
    var j = 0
    var temp = 0
    while (i < numElements-1){
      j = rand.nextInt(numElements-i-1) + i + 1
      temp = a(i)
      a(i) = a(j)
      a(j) = temp
      i += 1
    }
    a
  }

  /**
   * Generate a random permutation in parallel.
   * @param nThreads Number of threads to participate in the permutation operation.
   * @param threads Thread pool consisting of at least nThreads threads.
   * @param randSeed Random seed to use for the permutation.
   * @param nElements Number of elements to permute
   * @return Tuple of permutated ordering and its inverse: invOrder(i) = j iff order(j) = i
   */
  def parallelRandomPermutation(nThreads: Int, threads: ExecutorService, randSeed: Int, nElements: Int)
  : (Array[Int], Array[Int]) = {
    val masterRand : Random = new Random(randSeed)
    val rands : Array[Random] = (0 until nThreads).map(_ => new Random(masterRand.nextLong())).toArray
    val orderingStartPoint : Array[Int] = new Array[Int](nThreads)
    val ordering : Array[Int] = new Array[Int](nElements)
    val invOrder : Array[Int] = new Array[Int](nElements)
    var i = 0
    var hasChanged = true

    // First local permute
    val tasks_firstLocalPermute = (0 until nThreads).map(threadID => new Callable[Array[Int]]{
      override def call() = {
        val startIndex: Int =                                            math.floor( threadID     .toDouble / nThreads.toDouble * nElements.toDouble).toInt
        val endIndex  : Int = if (threadID == nThreads-1) nElements else math.floor((threadID + 1).toDouble / nThreads.toDouble * nElements.toDouble).toInt
        KnuthPermute((startIndex until endIndex).toArray, rands(threadID))
      }
    })
    val firstPermuteArrays : Array[Array[Int]] = threads.invokeAll(tasks_firstLocalPermute).map(_.get).toArray

    // Collect across processors, then second permute
    val tasks_collectAndPermute = (0 until nThreads).map(threadID => new Callable[Array[Int]]{
      override def call() = {
        var i = 0
        var j = 0
        var bucketPointer = 0
        var bucketNumElements = 0
        var localNumElements = 0
        var startIndex = 0
        var endIndex = 0
        while (i < nThreads){
          bucketPointer = (i + threadID) % nThreads
          bucketNumElements = firstPermuteArrays(bucketPointer).length
          startIndex =                                                    math.floor( threadID     .toDouble / nThreads.toDouble * bucketNumElements.toDouble).toInt
          endIndex   = if (threadID == nThreads-1) bucketNumElements else math.floor((threadID + 1).toDouble / nThreads.toDouble * bucketNumElements.toDouble).toInt
          localNumElements += (endIndex - startIndex)
          i += 1
        }
        val buffer : Array[Int] = new Array[Int](localNumElements)
        i = 0
        var counter = 0
        while (i < nThreads){
          bucketPointer = (i + threadID) % nThreads
          bucketNumElements = firstPermuteArrays(bucketPointer).length
          startIndex =                                                    math.floor( threadID     .toDouble / nThreads.toDouble * bucketNumElements.toDouble).toInt
          endIndex   = if (threadID == nThreads-1) bucketNumElements else math.floor((threadID + 1).toDouble / nThreads.toDouble * bucketNumElements.toDouble).toInt
          j = startIndex
          while (j < endIndex){
            buffer(counter) = firstPermuteArrays(bucketPointer)(j)
            j += 1
            counter += 1
          }
          i += 1
        }
        KnuthPermute(buffer.toArray, rands(threadID))
      }
    })
    val secondPermuteArrays : Array[Array[Int]] = threads.invokeAll(tasks_collectAndPermute).map(_.get).toArray

    // Create ordering
    i = 1
    while (i < nThreads){
      orderingStartPoint(i) = orderingStartPoint(i-1) + secondPermuteArrays(i-1).length
      i += 1
    }
    val tasks_ordering = (0 until nThreads).map(threadID => new Callable[Boolean]{
      override def call() = {
        var hasChanged = false
        var index = 0
        val numElements = secondPermuteArrays(threadID).length
        try {
          while (index < numElements) {
            if (ordering(orderingStartPoint(threadID) + index) != secondPermuteArrays(threadID)(index)) {
              ordering(orderingStartPoint(threadID) + index) = secondPermuteArrays(threadID)(index)
              hasChanged = true
            }
            index += 1
          }
          hasChanged
        }catch{
          case e: Exception =>{
            println(s"Thread $threadID, Error @ index = $index," +
              s"orderingStartPoint($threadID) = ${orderingStartPoint(threadID)}, ordering.length = ${ordering.length}, secondPermuteArrays($threadID).length = ${secondPermuteArrays(threadID).length}")
            e.printStackTrace()
            System.exit(-1)
            true
          }
        }
      }
    })
    hasChanged = true
    while (hasChanged) hasChanged = threads.invokeAll(tasks_ordering).map(_.get).reduce(_ | _)

    // Create invOrder
    val tasks_invOrder = (0 until nThreads).map(threadID => new Callable[Boolean]{
      override def call() = {
        var hasChanged = false
        var index = threadID
        while (index < nElements){
          if (invOrder(ordering(index)) != index){
            invOrder(ordering(index)) = index
            hasChanged = true
          }
          index += nThreads
        }
        hasChanged
      }
    })
    hasChanged = true
    while (hasChanged) hasChanged = threads.invokeAll(tasks_invOrder).map(_.get).reduce(_ | _)

    (ordering, invOrder)
  }

  /**
   * Computes the objective value of a clustering assignment to a given graph.
   * @param nThreads Number of threads to use for this check.
   * @param threads Thread pool consisting of at least nThreads threads.
   * @param graph Graph on which the correlation clustering has been done.
   * @param clusterID Assignment of vertices to clusters.
   * @param clusterID_atomic Atomic array corresponding to clusterID. This is taken as the ground truth.
   * @return correlation clustering objective value
   */
  def computeObjective(nThreads: Int,
                       threads: ExecutorService,
                       graph: SparseGraph,
                       clusterID : Array[Int],
                       clusterID_atomic: AtomicIntegerArray) : Long = {
    val nVertices: Int = graph.numNodes()

    val tasks = (0 until nThreads).map(threadID => new Callable[Long] {
      override def call() = {
        var objVal : Long = 0
        var index = threadID
        var v = 0
        var vClusterID = 0
        var ui = 0
        var u = 0
        var uClusterID = 0
        var uNSucc = 0
        var wi = 0
        var nSucc = 0
        val clusterHash : mutable.HashSet[Int] = new mutable.HashSet[Int]()
        while (index < nVertices){
          v = index
          if (clusterID(v)==0) clusterID(v) = clusterID_atomic.get(v)
          vClusterID = abs(clusterID(v))
          nSucc = graph.nSuccessor(v)
          // First, count bad +ve edges around this vertex
          ui = 0
          while (ui < nSucc){
            u = graph.succ(v, ui)
            if (clusterID(u)==0) clusterID(u) = clusterID_atomic.get(u)
            uClusterID = abs(clusterID(u))
            if (vClusterID != uClusterID){
              // Share edge but different cluster => bad +ve
              objVal += 1
            }
            ui += 1
          }
          // If center, count the bad -ve edges in the cluster
          if (clusterID(v) > 0){
            clusterHash.clear()
            clusterHash.add(v)
            ui = 0
            while (ui < nSucc){
              u = graph.succ(v, ui)
              if (clusterID(u)==0) clusterID(u) = clusterID_atomic.get(u)
              uClusterID = abs(clusterID(u))
              if (uClusterID == vClusterID) clusterHash.add(u)
              ui += 1
            }
            val clusterSize = clusterHash.size - 1
            ui = 0
            while (ui < nSucc){
              u = graph.succ(v, ui)
              if (clusterHash.contains(u)){
                wi = 0
                uNSucc = graph.nSuccessor(u)
                var clusterNeighbors = 0
                while (wi < uNSucc){
                  if (clusterHash.contains(graph.succ(u, wi))) clusterNeighbors += 1
                  wi += 1
                }
                objVal += (clusterSize - clusterNeighbors)
              }
              ui += 1
            }
          }
          index += nThreads
        }
        objVal
      }
    })
    threads.invokeAll(tasks).map(_.get()).reduce(_+_) / 2
  }

  /**
   * Given two sets of cluster IDs, check if they are equivalent
   * @param nThreads Number of threads to use for this check.
   * @param threads Thread pool consisting of at least nThreads threads
   * @param cid1 First set of cluster IDs.
   * @param cidA1 Atomic array corresponding to first set of cluster IDs. This is taken as the ground truth.
   * @param cid2 Second set of cluster IDs.
   * @param cidA2 Atomic array corresponding to the second set of cluster IDs. This is taken as the ground truth.
   * @return True iff the two sets of cluster IDs are the same.
   */
  def checkEquivalence(nThreads: Int,
                       threads: ExecutorService,
                       cid1: Array[Int], cidA1: AtomicIntegerArray,
                       cid2: Array[Int], cidA2: AtomicIntegerArray) : Boolean = {
    val nElements = cid1.length
    if (cid2.length != nElements || cidA1.length() != nElements || cidA2.length() != nElements){
      System.out.println(s"Lengths of cluster ID arrays do not match")
      false
    }else {
      val tasks_check = (0 until nThreads).map(threadID => new Callable[Boolean] {
        override def call() = {
          var isEquivalent = true
          var index = 0
          var id1 = 0
          var id2 = 0
          while (index < nElements){
            id1 = if (cid1(index)==0) cidA1.get(index) else cid1(index)
            id2 = if (cid2(index)==0) cidA2.get(index) else cid2(index)
            if (id1 != id2){
              System.out.println(s"cid1($index) == $id1 != $id2 == cid2($index)")
              isEquivalent = false
            }
            index += nThreads
          }
          isEquivalent
        }
      })
      threads.invokeAll(tasks_check).map(_.get).reduce(_|_)
    }
  }

}