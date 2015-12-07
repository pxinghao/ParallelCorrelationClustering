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

import java.util.concurrent.ExecutorService
import java.util.concurrent.atomic.{AtomicIntegerArray, AtomicInteger}

/**
 * SequentialRandomGreedyPeeling runs the serial KwikCluster algorithm of [10] Ailon, Charikar, and Newman (2008).
 * Given a graph and an ordering, KwikCluster sequentially processes the vertices in order.
 * When the algorithm encounters an unclustered vertex, it creates a new cluster, and puts all the unassigned neighbors of the new cluster center into that cluster.
 * @param graph Graph on which to run KwikCluster.
 */
class SequentialRandomGreedyPeeling(graph: SparseGraph) {
  /**
   * Actual function that executes KiwkCluster.
   * @param nThreads This parameter is ignored, since KwikCluster is a serial algorithm that is run on a single thread.
   * @param threads This parameter is ignored, since KwikCluster is a serial algorithm that is run on a single thread.
   * @param ordering The order in which to process the vertices: order(i) is the ith vertex to be processed
   * @param invOrder The inverse function of order: invOrder(v) is equal to order(i) iff order(i) = v
   * @param clusterID Output array: |clusterID(v)| is the cluster to which vertex v belongs; the id is positive iff v is a cluster center.
   * @return Number of edges touched, and number of BSP rounds. For serial KwikCluster, the number of BSP rounds is always set as 0.
   */
  def run(nThreads: Int, threads: ExecutorService, ordering: Array[Int], invOrder: Array[Int], clusterID: Array[Int])
  : (Long, Int) = {
    val nElements : Int = graph.numNodes()

    var numEdgesExamined : Long = 0

    var index = 0
    var v = 0
    var u = 0
    var ui = 0
    var nSucc = 0
    while (index < nElements){
      v = ordering(index)
      if (clusterID(v) == 0){
        // New cluster
        clusterID(v) = index+1
        // Peel
        ui = 0
        nSucc = graph.nSuccessor(v)
        while (ui < nSucc){
          u = graph.succ(v, ui)
          if (clusterID(u) == 0){
            clusterID(u) = -(index+1)
          }
          ui += 1
        }
        numEdgesExamined += nSucc
      }
      index += 1
    }
    (numEdgesExamined, 0)
  }


}
