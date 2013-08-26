package abra.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

public class CalcMismatchDensity {

	private String ref;
	private String baseDir;
	private String inputFile;
	private String currPid = "";
	private String origBam;
	private String abraBam;
	private CompareToReference3 c2r;
	private int minBaseQuality = 20;
	private int readLength;
	
	public CalcMismatchDensity(String ref, String baseDir, String inputFile, String origBam, String abraBam, int readLength) {
		this.ref = ref;
		this.baseDir = baseDir;
		this.inputFile = inputFile;
		this.origBam = origBam;
		this.abraBam = abraBam;
		this.readLength = readLength;
	}
	
	public void run() throws IOException {
		
		System.err.println("Loading reference");
		c2r = new CompareToReference3();
		c2r.init(ref);
		c2r.setMinBaseQuality(20);
		
		BufferedReader reader = new BufferedReader(new FileReader(inputFile));
		
		SAMFileReader abraReader = null;
		SAMFileReader origReader = null;
		
		long s = System.currentTimeMillis();
		
		String line = reader.readLine();
		while (line != null) {
			String[] fields = line.split("\t");
			String pid = fields[0];
			if (!pid.equals(currPid)) {
				currPid = pid;
				System.err.println(currPid);
				abraReader = initReader(abraReader, baseDir + "/" + pid + "/" + abraBam);
				origReader = initReader(origReader, baseDir + "/" + pid + "/" + origBam);
			}
			
			String chromosome = fields[1];
			int pos = Integer.parseInt(fields[2]);
			int refLen = Integer.parseInt(fields[4]);
			int altLen = Integer.parseInt(fields[5]);
			int length = 0;
			CigarOperator indelType;
			
			if (refLen > 1) {
				indelType = CigarOperator.D;
				length = refLen - 1;
			} else if (altLen > 1) {
				indelType = CigarOperator.I;
				length = altLen - 1;
			} else {
				throw new IllegalArgumentException("Not an indel: " + line);
			}
			
			double abraMd = calcMismatchDensity(abraReader, indelType, chromosome, pos, length);
			double origMd = calcMismatchDensity(origReader, indelType, chromosome, pos, length);
			
			System.out.println(line + "\t" + abraMd + "\t" + origMd);
			
			line = reader.readLine();
		}
		
		reader.close();
		
		long e = System.currentTimeMillis();
		
//		System.out.println("Elapsed: " + (e-s));
	}
	
	private double calcMismatchDensity(SAMFileReader reader, CigarOperator indelType, String chromosome, int pos, int length) {
		
		int numMismatches = 0;
		int numBases = 0;
		
		CloseableIterator<SAMRecord> iter =  reader.queryOverlapping(chromosome, pos-readLength, pos+length+readLength);

		while (iter.hasNext()) {
			SAMRecord read = (SAMRecord) iter.next();
			
			String yo = (String) read.getAttribute("YO");

			if ((yo == null) || (!yo.equals("N/A"))) {
				numMismatches += c2r.noiseAroundIndel(read, indelType, pos, length);
				
				for (int qual : read.getBaseQualities()) {
	//				qual = qual - '!';
					if (qual >= minBaseQuality) {
						numBases += 1;
					}
				}
			}
		}

		iter.close();
		
		return numBases == 0 ? 0.0 : (double) numMismatches / (double) numBases;
	}
	
	private SAMFileReader initReader(SAMFileReader old, String input) {
		if (old != null) {
			old.close();
		}
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);

		return reader;
	}
	
	public static void main(String[] args) throws Exception {
		String ref = args[0];
		String baseDir = args[1];
		String inputFile = args[2];
		String origBam = args[3];
		String abraBam = args[4];
		int readLength = Integer.parseInt(args[5]);

//		String ref = "/home/lmose/reference/chr1/1.fa";
//		String baseDir = "/home/lmose/dev/ayc/mismatch_density";
//		String inputFile = "/home/lmose/dev/ayc/mismatch_density/gl.txt";
//		String origBam = "n1.bam";
//		String abraBam = "n.abra.1.bam";
//		int readLength = 100;
		
		CalcMismatchDensity cmd = new CalcMismatchDensity(ref, baseDir, inputFile, origBam, abraBam, readLength);
		
		cmd.run();
		
		System.err.println("Done.");
	}
}
