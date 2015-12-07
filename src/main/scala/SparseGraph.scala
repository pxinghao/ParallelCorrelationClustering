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

import it.unimi.dsi.webgraph.BVGraph

/**
  * SparseGraph representation
  */
class SparseGraph() {
  /**
   * Number of vertices
   */
  private var nVertices = 0

  /**
   * Number of edges
   */
  private var nEdges = 0

  /**
   * Offset of each vertex into edge list
   */
  private var rowOffsets = new Array[Int](0)

  /**
   * Number of successors / neighbors of each vertex
   */
  private var rowNSucc   = new Array[Int](0)

  /**
   * Edge list
   */
  private var edgeList   = new Array[Int](0)

  /**
   * Populates the internal data structures.
   * @param m_nVertices  Number of vertices
   * @param m_nEdges     Number of edges
   * @param m_rowOffsets Offset of each vertex into edge list
   * @param m_rowNSucc   Number of successors / neighbors of each vertex
   * @param m_edgeList   Edge list
   */
  def create(m_nVertices: Int, m_nEdges: Int, m_rowOffsets: Array[Int], m_rowNSucc: Array[Int], m_edgeList: Array[Int]) = {
    nVertices  = m_nVertices
    nEdges     = m_nEdges
    rowOffsets = m_rowOffsets
    rowNSucc   = m_rowNSucc
    edgeList  = m_edgeList
  }

  /**
   * Load a BVGraph into SparseGraph format
   * @param graph BVGraph to load
   */
  def loadFromBVGraph(graph: BVGraph) : Unit = {
    nVertices = graph.numNodes()
    nEdges = 0
    rowOffsets = new Array[Int](nVertices)
    rowNSucc   = new Array[Int](nVertices)
    edgeList  = new Array[Int](graph.numArcs().toInt)

    var i = 0
    while (i < nVertices){
      rowOffsets(i) = nEdges
      var nSucc = 0
      val siter = graph.successors(i)
      var j = siter.nextInt()
      while (j != -1){
        if (j != i){
          edgeList(nEdges) = j
          nEdges += 1
          nSucc += 1
        }
        j = siter.nextInt()
      }
      rowNSucc(i) = nSucc
      i += 1
    }

  }

  /**
   * Number of successors / neighbors of a vertex
   * @param i Vetex index
   * @return  Number of successors
   */
  def nSuccessor(i: Int) = {
    rowNSucc(i)
  }

  /**
   * Retrieves a particular neighbor of a given vetex
   * @param i Vertex queried
   * @param j Index of neighbor to retrieve
   * @return  jth neighbor of ith vertex
   */
  def succ(i: Int, j: Int) = {
    edgeList(rowOffsets(i) + j)
  }

  /**
   * Number of vertices
   * @return Number of vertices
   */
  def numNodes() = {
    nVertices
  }

  /**
   * Number of edges
   * @return Number of edges
   */
  def numArcs() = {
    nEdges
  }


}
