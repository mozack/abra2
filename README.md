# ABRA - Assembly Based ReAligner

## Introduction

ABRA is a realigner for next generation sequencing data.  It uses localized assembly and global realignment to align reads more accurately, thus improving downstream analysis (detection of indels and complex variants in particular).

Here is an ABRA realigned region (original reads on top, ABRA realigned reads on bottom):

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

The above command allocates 4GB for the java heap.  ABRA includes native code that runs outside of the JVM.  We have found that for typical exome processing using 8 threads, 16GB is more than sufficient.


### Parameters
parameter | value
------ | -------
--in | One or more input BAMs delimited by comma
--out | One or more output BAM's corresponding to the set of input BAMs
--ref  | BWA indexed reference genome.
--kmer | Comma delimited list of kmers used for assembly.  Smallest value is used by default, larger values are used if necessary.
--targets | BED file describing target assembly regions (Usually corresponds to capture targets)
--working | Temp working directory

### Somatic  mode

If working with tumor/normal pairs, it is highly recommended to assemble your samples together.  To do this, simply specify multiple input and output BAM files on the command line. 

### Output
ABRA produces one or more realigned BAMs.  It is currently necessary to sort and index the output.  At present, the mate information may not be 100% accurate.  Samtools fixmate or Picard Tools FixMateInformation may be used to correct this.

Reads that have been realigned will contain a YO tag indicating their original alignment position.  Reads that were originally unaligned will have a YO value of N/A.

After the output BAM file is sorted and indexed it can be passed into a variant caller such as [FreeBayes](https://github.com/ekg/freebayes)

### Demo / test data
A test data set and example command line is available under the demo directory.  Edit demo.bash to specify your hg19 reference location (indexed by a recent version of bwa) and run demo.bash.