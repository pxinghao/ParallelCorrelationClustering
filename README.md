#Parallel Correlation Clustering

This repository contains the code for the parallel algorithms presented in our <a href="paper/pan-etal.nips2015corrclus1.pdf">NIPS 2015 paper</a> for parallel correlation clustering.

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


<H3>Running examples using SBT</H3>

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
