abra
====

Assembly Based ReAligner

Running ABRA currently requires bwa 0.7.5a in the command path

Sample command line for v0.69:

```
java -Xss8M -Xmx20G -XX:MaxPermSize=256M  -jar $JAR --in $BAM --kmer 43,53,63,73,83 --mc-mapq 25 --mcl 102 --mcr -1.0 --mnf 2 --umnf 2 --mpc 50000 --out $ABRA_BAM --ref $REF --targets wxs.gtf --threads 8 --working $WORK --mur 50000000 --paired --no-unalign --mbq 5
```

Pre-package jars are available for 64 bit linux: https://github.com/mozack/abra/releases

The jar does contain native code, so it may be necessary to re-compile for your platform.  We have tested on:
64 bit CentOS 5.8
64 bit CentOS 6.2
64 bit Red Hat 5.5

Building ABRA requires JDK 7, Maven and g++
Just run make

More details to come...

Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved.
