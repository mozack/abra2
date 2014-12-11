# ABRA - Assembly Based ReAligner

## Introduction

ABRA is a realigner for next generation sequencing data.  It uses localized assembly and global realignment to align reads more accurately, thus improving downstream analysis (detection of indels and complex variants in particular).

Here is an ABRA realigned region (original reads on top, ABRA realigned reads on bottom).  The original set of reads have rather "noisy" alignments with several variations from the reference and a fair bit of high quality soft clipping.  The ABRA realignments present a more parsimonious representation of the reads including a previously unobserved large deletion. 

![ABRA Example](https://raw.githubusercontent.com/mozack/abra/master/misc/example.png)

## Building

Building ABRA requires JDK 7, Maven and g++

Just run make.  An executable jar will be generated under the target directory.  i.e. target/abra-0.77-SNAPSHOT-jar-with-dependencies.jar

Pre-built jars are available for 64 bit linux here: https://github.com/mozack/abra/releases

Note, that the jar does contain native code.  While we have tested on a variety of platforms, we cannot guarantee that the pre-built jar will work everywhere.  In some cases, it may be necessary to compile for your specific platform.

## Running

Running ABRA currently requires bwa 0.7.9a (or similar) in the command path and a recent version of Java.

Sample command line for v0.86:

```
java -Xmx4G -jar $JAR --in input.bam --out output.bam --ref hg19.fasta --targets targets.bed --threads 8 --working abra_temp_dir > abra.log 2>&1
```

The above command allocates 4GB for the java heap.  ABRA includes native code that runs outside of the JVM.  We have found that for typical exome processing using 8 threads, 16GB total is more than sufficient.

Note: The code at the HEAD may occasionally be unstable.  It is recommended to work from a release.

### Parameters
parameter | value
------ | -------
--in | One or more input BAMs delimited by comma
--out | One or more output BAM's corresponding to the set of input BAMs
--ref  | BWA indexed reference genome.
--targets | BED file describing target assembly regions (Usually corresponds to capture targets) with optional kmer lengths for each region.  Targets must be sorted by position in ascending order within each chromosome.
--working | Temp working directory

### Somatic  mode

If working with tumor/normal pairs, it is highly recommended to assemble your samples together.  To do this, simply specify multiple input and output BAM files on the command line.

In somatic mode, you may also consider using option ```--lr repeat_file```  This option allows for detection of moderate length repeats that could not be resolved at nucleotide precision.  This feature is experimental.

Additionally, translocations may be detected by using ```-sv translocation_file```  The output file will contain putative breakpoints with the number of reads aligned to each breakpoint in the output file.  This feature is experimental.

Sample usage:
```
java -Xmx4G -jar $JAR --in normal.bam,tumor.bam --out normal.abra.bam,tumor.abra.bam --ref hg19.fasta --targets targets.bed --threads 8 --working abra_temp_dir > abra.log 2>&1
```

### Output
ABRA produces one or more realigned BAMs.  It is currently necessary to sort and index the output.  At present, the mate information may not be 100% accurate.  Samtools fixmate or Picard Tools FixMateInformation may optionally be used to correct this.

Reads that have been realigned will contain a YO tag indicating their original alignment position.  Reads that were originally unaligned will have a YO value of N/A.

After the output BAM file is sorted and indexed it can be passed into a variant caller such as [FreeBayes](https://github.com/ekg/freebayes) for germline calling and [Cadabra](https://github.com/mozack/abra#cadabra), [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller) or [UNCeqr](http://lbg.med.unc.edu/~mwilkers/unceqr_dist) for somatic calling.

### Demo / test data
A test data set and example command line is available under the demo directory.  Edit demo.bash to specify your hg19 reference location (indexed by a recent version of bwa) and run demo.bash.

### Cadabra

Cadabra is a somatic indel caller that works specifically with ABRA alignments.  It uses SAM tags embedded by ABRA to identify assembled indels with support in the tumor, and a general lack of support in the normal.  It performs very well according to our tests.

Quality scores are calculated using a Phred scaled Fisher's Exact Test.

Cadabra is embedded within the ABRA jar file.

Sample usage:
```
java -Xmx8G -cp $JAR abra.cadabra.Cadabra hg19.fasta normal.abra.bam tumor.abra.bam > cadabra.vcf
```

