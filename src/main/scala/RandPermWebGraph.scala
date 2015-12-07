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

import java.io.{File, PrintWriter}
import java.util.Random

import it.unimi.dsi.webgraph.BVGraph

/**
 * Randomly permutes a BVGraph
 */
object RandPermWebgraph {
  def main(args: Array[String]) = {
    val argmap = args.map { a =>
      val argPair = a.split("=")
      val name = argPair(0).toLowerCase
      val value = argPair(1)
      (name, value)
    }.toMap

    // Default options
    val basename        = argmap.getOrElse("inputfile",  "data/eswiki-2013")
    val randSeed        = argmap.getOrElse("randseed",   "43").toInt

    val gen = new Random(randSeed)

    val g = BVGraph.load(basename)
    val nVertices = g.numNodes()

    // randPerm maps old to new : randPerm(old) = new
    // revIndex maps new to old : revIndex(new) = old

    val randPerm = (0 until nVertices).toArray
    shuffle(randPerm, gen)

    val revIndex = new Array[Int](nVertices)
    var ii = 0
    while (ii < nVertices){
      //println(s"$ii -> ${randPerm(ii)}")
      revIndex(randPerm(ii)) = ii
      ii += 1
    }

    val asciiWriter = new PrintWriter(new File(basename + "_rand.graph-txt"))
    asciiWriter.print(s"$nVertices\n")
    var i = 0
    while (i < nVertices){
      var adj = g.successorArray(revIndex(i))
      var j = 0
      while (j < adj.length){
        adj(j) = randPerm(adj(j))
        j += 1
      }
      adj = adj.sorted
      j = 0
      while (j < adj.length - 1){
        asciiWriter.print(adj(j) + " ")
        j += 1
      }
      if (adj.length > 0) asciiWriter.print(adj(j))
      asciiWriter.println()
      i += 1
    }
    asciiWriter.close()

    BVGraph.main(Array("-g", "ASCIIGraph", basename + "_rand", basename + "_rand"))

    BVGraph.main(Array("-o", "-O", "-L", basename + "_rand"))
  }

  def shuffle[T](array: Array[T], rnd: Random): Array[T] = {
    for (n <- Iterator.range(array.length - 1, 0, -1)) {
      val k = rnd.nextInt(n + 1)
      val t = array(k); array(k) = array(n); array(n) = t
    }
    array
  }
}
