# Path to ABRA jar file (with dependencies)
# If you've downloaded the jar, set this to the appropriate location
# Compiled jars will be under the target directory
JAR=../target/abra-0.75-SNAPSHOT-jar-with-dependencies.jar

# Path to hg19 reference (indexed by BWA.  bwa mem must be in your command path)
#REF=/datastore/rclbg/nextgenout3/MOSE_TEST/abra/brca/ref/GRCh37-lite.fa
REF=<path to hg19 reference>

java -Xmx4G -jar $JAR --ref $REF --in abra.demo.bam --out abra.demo.realigned.bam --kmer 43,53,63,73,83 --working abra_temp_dir --targets demo.bed