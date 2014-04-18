abra
====

Assembly Based ReAligner

Running ABRA currently requires bwa 0.7.5a (or similar) in the command path

Sample command line for v0.75:

```
java -Xmx4G -jar $JAR --in normal.bam,tumor.bam --kmer 43,53,63,73,83 --out normal.abra.bam,tumor.abra.bam --ref $REF --targets wxs.bed --threads 8 --working abra.work > abra.log 2>&1
```

Pre-packaged jars are available for 64 bit linux here: https://github.com/mozack/abra/releases

The jar does contain native code.  While we have tested on a variety of platforms, it may be necessary to re-compile for your own.

Building ABRA requires JDK 7, Maven and g++<br/>
Just run make.  An executable jar will be generated under the target directory.  i.e. target/abra-0.75-SNAPSHOT-jar-with-dependencies.jar

More details to come....

Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved.
