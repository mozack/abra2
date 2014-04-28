# ABRA - Assembly Based ReAligner

## Introduction

ABRA is a realigner for next generation sequencing data.  It uses localized assembly and global realignment to align reads more accurately, thus improving downstream analysis (detection of indels and complex variants in particular).

Here is an example of ABRA realigned reads (original reads on top, ABRA realigned reads on bottom):

![ABRA Example](https://raw.githubusercontent.com/mozack/abra/master/misc/example.png)

## Building

Building ABRA requires JDK 7, Maven and g++

Just run make.  An executable jar will be generated under the target directory.  i.e. target/abra-0.75-SNAPSHOT-jar-with-dependencies.jar

Pre-built jars are available for 64 bit linux here: https://github.com/mozack/abra/releases

Note, that the jar does contain native code.  While we have tested on a variety of platforms, we cannot guarantee it will work everywhere.

## Running

Running ABRA currently requires bwa 0.7.5a (or similar) in the command path and a recent version of Java.

Sample command line for v0.75:

```
java -Xmx4G -jar $JAR --in input.bam --kmer 43,53,63,73,83 --out output.bam --ref hg19.fasta --targets targets.bed --threads 8 --working abra_temp_dir > abra.log 2>&1
```

### Parameters
--in <input BAM(s)>	| (Multiple BAMs may be specified, separated by comma)
--out <output BAM(s)> | (The number of output BAMs must match the number of input BAMs)
--ref <reference fasta> | (Must be indexed by BWA)
--kmer <comma delimited list of kmers for assembly> | (Smallest kmer is used by default, larger kmers are used when necessary)
--targets <BED file describing target assembly regions> | (Usually corresponds to capture targets)
--working | <temp directory>

### Somatic  mode

If working with tumor/normal pairs, it is highly recommended to assemble your samples together.  To do this, simply specify multiple input and output BAM files on the command line. 

### Output
ABRA produces 1 or more realigned BAMs.  It is currently necessary to sort and index the output.  At present, the mate information may not bee 100% accurate.  Samtools fixmate or Picard Tools FixMateInformation may be used to correct this.

Reads that have been realigned will contain a YO tag indicating their original alignment position.  Reads that were originally unaligned will have a YO value of N/A.

After the output BAM file is sorted and indexed it can be passed into a variant caller such as [FreeBayes](https://github.com/ekg/freebayes)