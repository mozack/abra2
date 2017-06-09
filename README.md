# ABRA2

ABRA2 is an updated implementation of [ABRA](https://github.com/mozack/abra) featuring:
* RNA support
* Improved scalability (Human whole genomes now supported)
* Improved accuracy
* Improved stability and usability (BWA is no longer required to run ABRA although we do recommend BWA as the initial aligner for DNA)

## Running

ABRA2 requires Java 8.

### DNA

Sample command for DNA:

```java -Xmx16G -jar abra2-2.07.jar --in normal.bam,tumor.bam --out normal.abra.bam,tumor.abra.bam --ref hg38.fa --threads 8 --targets targets.bed --dist 1000 --tmpdir /your/tmpdir > abra.log```

The above accepts normal.bam and tumor.bam as input and outputs sorted realigned BAM files named normal.abra.bam and tumor.abra.bam

* Input files must be sorted by coordinate and index
* Output files are sorted
* The tmpdir may grow large.  Be sure you have sufficient space there (at least equal to the input file size)
* The targets argument is not required.  When omitted, the entire genome will be eligible for realignment.

### RNA

Abra2 is capable of utilizing junction information to aid in assembly and realignment.  It has been tested only on STAR output to date.

Sample command for RNA:

```java -Xmx16G -jar abra2-2.07.jar --in star.bam --out star.abra.bam --ref hg38.fa --junctions SJ.out.tab --threads 8 --gtf gencode.v26.annotation.gtf --dist 500000 --tmpdir /your/tmpdir  > abra2.log 2>&1```

Here, star.bam is the input bam file and star.abra.bam is the output bam file.

Junctions observed during alignment can be passed in using the ```--junctions``` param.  This corresponds to the SJ.out.tab file output by STAR.

Annotated junctions can be passed in using the ```--gtf``` param.  See: https://www.gencodegenes.org/releases/current.html  
It is beneficial to use both of the junction related options.

The software is currently considered beta quality.

