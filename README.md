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

Running ABRA currently requires bwa 0.7.5a (or similar) in the command path and a recent version of Java.

Sample command line for v0.77:

```
java -Xmx4G -jar $JAR --in input.bam --kmer 43,53,63,73,83 --out output.bam --ref hg19.fasta --targets targets.bed --threads 8 --working abra_temp_dir > abra.log 2>&1
```

The above command allocates 4GB for the java heap.  ABRA includes native code that runs outside of the JVM.  We have found that for typical exome processing using 8 threads, 16GB total is more than sufficient.

Note: The code at the HEAD may occasionally be unstable.  It is recommended to work from a release.

### Parameters
parameter | value
------ | -------
--in | One or more input BAMs delimited by comma
--out | One or more output BAM's corresponding to the set of input BAMs
--ref  | BWA indexed reference genome.
--kmer | Comma delimited list of kmer lengths used for assembly.  Smallest value is used by default, larger values are used if necessary.  Ignored if kmer sizes are supplied in the target bed file.
--targets | BED file describing target assembly regions (Usually corresponds to capture targets) with optional kmer lengths for each region
--working | Temp working directory

### Kmer length selection

As of version 0.77, kmer sizes can by determined based upon the reference content using KmerSizeEvaluator.  We have seen improved results using this method.

Usage:
```
java -Xmx4G -cp abra.jar abra.KmerSizeEvaluator <read_length> <reference_fasta> <output_bed> <num_threads> <input_bed> <working_dir>
```

This will create file "output_bed" which contains the regions specified by "input_bed" with an additional column for kmer size to be used for that region.
The "output_bed" file can then be passed as input to ABRA via the --targets option.

Example:
```
java -Xmx4G -cp abra.jar abra.KmerSizeEvaluator 100 ref/ucsc.hg19.fasta abra_kmers.bed 8 wxs.bed abra_kmer_temp > kmers.log 2>&1
```

### Somatic  mode

If working with tumor/normal pairs, it is highly recommended to assemble your samples together.  To do this, simply specify multiple input and output BAM files on the command line.

In somatic mode, you may also consider using option ```--lr repeat_file```  This option allows for detection of moderate length repeats that could not be resolved at nucleotide precision.  This feature is experimental. 

### Output
ABRA produces one or more realigned BAMs.  It is currently necessary to sort and index the output.  At present, the mate information may not be 100% accurate.  Samtools fixmate or Picard Tools FixMateInformation may optionally be used to correct this.

Reads that have been realigned will contain a YO tag indicating their original alignment position.  Reads that were originally unaligned will have a YO value of N/A.

After the output BAM file is sorted and indexed it can be passed into a variant caller such as [FreeBayes](https://github.com/ekg/freebayes), [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller) and [UNCeqr](http://lbg.med.unc.edu/~mwilkers/unceqr_dist)

### Demo / test data
A test data set and example command line is available under the demo directory.  Edit demo.bash to specify your hg19 reference location (indexed by a recent version of bwa) and run demo.bash.

