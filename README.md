# ABRA2

ABRA2 is an updated implementation of ABRA (https://github.com/mozack/abra) featuring:
* RNA support
* Improved scalability (Human whole genomes now supported)
* Improved accuracy
* Improved stability and usability (BWA is no longer required to run ABRA although we do recommend BWA as the initial aligner for DNA)

## Running

ABRA2 requires Java 8.

### RNA

Abra2 is capable of utilizing junction information to aid in assembly and realignment.  It has been tested only on STAR output to date.

Sample command for RNA:

```java -Xmx16G -jar abra2-0.01.jar --in star.bam --out star.abra.bam --ref hg19.fa --junctions SJ.out.tab --threads 8 --gtf gencode.v19.annotation.gtf  > abra2.log 2>&1```

Here, star.bam is the input bam file and star.abra.bam is the output bam file.

Junctions observed during alignment can be passed in using the ```--junctions``` param.  This corresponds to the .tab file output by STAR.

Annotated junctions can be passed in using the ```--gtf``` param.  It is beneficial to use both of the junction related options.

The software is currently considered beta quality. 
