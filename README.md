# ABRA2

ABRA2 is an updated implementation of [ABRA](https://github.com/mozack/abra) featuring:
* RNA support
* Improved scalability (Human whole genomes now supported)
* Improved accuracy
* Improved stability and usability (BWA is no longer required to run ABRA although we do recommend BWA as the initial aligner for DNA)

Manuscript: https://doi.org/10.1093/bioinformatics/btz033

## Running

ABRA2 requires Java 8.
We recommend running from a pre-compiled release.
Go to the Releases tab to download a recent version.

### DNA

Sample command for DNA:

```java -Xmx16G -jar abra2.jar --in normal.bam,tumor.bam --out normal.abra.bam,tumor.abra.bam --ref hg38.fa --threads 8 --targets targets.bed --tmpdir /your/tmpdir > abra.log```

The above accepts normal.bam and tumor.bam as input and outputs sorted realigned BAM files named normal.abra.bam and tumor.abra.bam

* Input files must be sorted by coordinate and index
* Output files are sorted
* The tmpdir may grow large.  Be sure you have sufficient space there (at least equal to the input file size)
* The targets argument is not required.  When omitted, the entire genome will be eligible for realignment.

### RNA

ABRA2 is capable of utilizing junction information to aid in assembly and realignment.  It has been tested only on STAR output to date.

Sample command for RNA:

```java -Xmx16G -jar abra2.jar --in star.bam --out star.abra.bam --ref hg38.fa --junctions bam --threads 8 --gtf gencode.v26.annotation.gtf --dist 500000 --sua --tmpdir /your/tmpdir  > abra2.log 2>&1```

Here, star.bam is the input bam file and star.abra.bam is the output bam file.

Junctions observed during alignment can be passed in using the ```--junctions``` param.  The input file format is similar to the SJ.out.tab file output by STAR.  If ```bam``` is specified, ABRA2 will dynamically identify splice junctions from the BAM file on the fly.  Note that the SJ.out.tab file contains only junctions deemed "high quality" by STAR.  The complete set of all splice junctions can be identified using the program ```abra.cadabra.SpliceJunctionCounter```

Annotated junctions can be passed in using the ```--gtf``` param.  See: https://www.gencodegenes.org/releases/current.html  
It is beneficial to use both of the junction related options.

Known indels can be passed in using the --in-vcf argument.  Unannotated junctions originally identified as splices by the aligner may be converted to deletions if a known deletion is matched.  Consider this option if you have indels detected from DNA for the same sample / subject.  It is not recommended to use large datasets when using this option (i.e. don't pass in dbSNP).


