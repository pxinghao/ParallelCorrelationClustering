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
 * Given a BVGraph, creates the offsets file, then randomly permutes and symmetrizes the graph
 */
object PreprocessGraph{
  def main(args: Array[String]) = {
    val argmap = args.map { a =>
      val argPair = a.split("=")
      val name = argPair(0).toLowerCase
      val value = argPair(1)
      (name, value)
    }.toMap

    // Default options
    val basename = argmap.getOrElse("inputfile",  "data/eswiki-2013")

    println("Creating offsets for graph.. ")
    BVGraph.main(Array("-o", "-O", "-L", basename))

    println("Randomly permute graph")
    RandPermWebgraph.main(Array("inputfile=" + basename))

    println("Symmetrize graph")
    SymmetrizeWebgraph.main(Array("inputfile=" + basename + "_rand"))

  }
}

