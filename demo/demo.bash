# Demonstration command line for ABRA.
# a recent version of java must be in your command path
# Assumes your current working directory is the demo directory

# Path to ABRA jar file (with dependencies)
# If you've downloaded the jar, set this to the appropriate location
# Compiled jars will be under the target directory
JAR=../target/abra2-*-jar-with-dependencies.jar

# Path to hg19 reference
#REF=/datastore/rclbg/nextgenout3/MOSE_TEST/abra/brca/ref/GRCh37-lite.fa
REF=<path to hg19 reference>

echo "ABRA demo starting..."

java -Xmx4G -jar $JAR --ref $REF --in abra_demo.bam --out abra_demo_realigned.bam --targets demo.bed > abra_demo.log 2>&1

echo "ABRA demo done.  Realigned BAM: abra_demo_realigned.bam."