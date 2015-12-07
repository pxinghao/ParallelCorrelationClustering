name := "ParallelCorrelationClustering"

version := "1.0"

scalaVersion := "2.11.6"


libraryDependencies ++= Seq(
  "org.apache.commons" % "commons-math3" % "3.2",
  "it.unimi.dsi" % "webgraph" % "3.4.0",
  "org.scalanlp" %% "breeze" % "0.10",
  "org.scalanlp" %% "breeze-natives" % "0.10"
)

resolvers ++= Seq(
  "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/snapshots/",
  "Sonatype Releases" at "https://oss.sonatype.org/content/repositories/releases/"
)


