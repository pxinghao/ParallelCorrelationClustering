#Parallel Correlation Clustering

This repository contains the code for the parallel algorithms presented in our <a href="paper/pan-etal.nips2015corrclus1.pdf">NIPS 2015 paper</a> for parallel correlation clustering.
We take a database-transactional view of the serial KwikCluster algorithm which provides an expected objective value of 3 OPT.
Each set of operations on a vertex is treated as a transaction to be executed in parallel with other transactions.

Our first algorithm, C4 (correlation clustering using concurrency control), resolves conflicts between transactions, but either forcing later transactions to wait, or correcting the clustering assignment of the later conflicting transaction.
As a result, we guarantee serial equivalence with KwikCluster, and specifically maintains the 3 OPT approximation.

The second algorithm, ClusterWild!, ignores potential conflicts between transactions and simply executes them in parallel.
Despite the conflicts, we are able to analytically show that ClusterWild! (BSP) obtains a reasonable error.
In practice, the objective value of ClusterWild! (BSP) is almost the same as that of KwikCluster.

Our code implements the BSP and asynchronous versions of C4 and ClusterWild!, in addition to the serial KwikCluster algorithm and a distributed algorithm (CDK) described in this <a href="http://dl.acm.org/citation.cfm?id=2623743&dl=ACM&coll=DL&CFID=540027751&CFTOKEN=13236671">paper</a>.

##Running examples without compilation
As a demonstration, <a href=src/main/scala/Example.scala>Example.scala</a> runs the various correlation clustering algorithms.
Data for the examples are provided in the <a href=data/>data</a> directory, sourced from <a href="http://law.di.unimi.it/datasets.php">Laboratory for Web Algorithms</a>.


To run Example,
<pre>java -cp bin/parrCorrClus.jar Example inputfile=&lt;filename&gt; maxNThreads=&lt;maximum number of threads&gt;</pre>
For example,
<pre>java -cp bin/parrCorrClus.jar Examples inputfile=data/eswiki-2013_rand_symm maxNThreads=4</pre>
The default options runs serial KwikCluster (using 1 thread), C4 and ClusterWild! (BSP and asynchronous versions, using 1 to maxNThreads threads), and the CDK algorithm (using maxNThreads).
For more options, see the documentation on <a href=api_doc/index.html#Example$>Example.scala</a>

<H3>Preprocessing graph data for examples</H3>
Graph data downloaded from <a href="http://law.di.unimi.it/datasets.php">Laboratory for Web Algorithms</a> can be pre-processed (randomizing order and symmetrizing graph) for use by running <a href=src/main/scala/Preprocess.scala>Preprocess.scala</a> with the command line argument
<pre>java -cp bin/parrCorrClus.jar PreprocessGraph inputfile=&lt;filename&gt;</pre>
For example,
<pre>java -cp bin/parrCorrClus.jar PreprocessGraph inputfile=data/eswiki-2013</pre>


##Running examples using SBT

If you want to modify the code and run the examples again you can either recompile using bin/update_bin or using the sbt launch scripts as below:

To run Example (this will recompile the code if necessary),
<pre>sbt/sbt "run-main Examples inputfile=&lt;filename&gt; maxNThreads=&lt;maximum number of threads&gt;"</pre>
For example,
<pre>sbt/sbt "run-main Examples inputfile=data/eswiki-2013_rand_symm maxNThreads=4"</pre>
For more options, see the documentation on <a href=api_doc/index.html#Example$>Example.scala</a>


<H3>Preprocessing graph data for examples</H3>
Graph data downloaded from <a href="http://law.di.unimi.it/datasets.php">Laboratory for Web Algorithms</a> can be pre-processed (randomizing order and symmetrizing graph) for use by running <a href=src/main/scala/Preprocess.scala>Preprocess.scala</a> with the command line argument
<pre>sbt/sbt "run-main PreprocessGraph inputfile=&lt;filename&gt;"</pre>
For example,
<pre>sbt/sbt "run-main PreprocessGraph inputfile=data/eswiki-2013"</pre>
