package abra;

import java.io.File;
import java.io.IOException;


import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Converts SAM/BAM file to FASTQ
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Sam2Fastq {
	
	private FastqOutputFile output1;
	private FastqOutputFile output2;
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	private boolean shouldIdentifyEndByReadId = false;
	private boolean isMapspliceFusions = false;
	private String end1Suffix;
	private String end2Suffix;

	/**
	 * Convert the input paired end SAM/BAM file into 2 fastq files.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 */
	public void convert(String inputSam, String outputFastq1, String outputFastq2) throws IOException {
		String last1Read = "";
		String last2Read = "";
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);

        output1 = new FastqOutputFile();
        output1.init(outputFastq1);
        output2 = new FastqOutputFile();
        output2.init(outputFastq2);
        
        int output1Count = 0;
        int output2Count = 0;
        int lineCnt = 0;
        
        for (SAMRecord read : reader) {
        	if (isFirstInPair(read)) {
        		if (!read.getReadName().equals(last1Read)) {
        			output1.write(samReadToFastqRecord(read));
        			last1Read = read.getReadName();
        			output1Count += 1;
        		}
        	} else if (isSecondInPair(read)) {
        		if (!read.getReadName().equals(last2Read)) {
        			output2.write(samReadToFastqRecord(read));
        			last2Read = read.getReadName();
        			output2Count += 1;
        		}
        	} else {
        		System.out.println("Don't know what to do with unpaired read: " + read.getSAMString());
        	}
        	
            lineCnt++;
            if ((lineCnt % 1000000) == 0) {
                System.out.println("record: " + lineCnt);
            }
        }
                
        output1.close();
        output2.close();
        reader.close();
        
        if (output1Count != output2Count) {
        	throw new IllegalStateException("Non-symmetrical read counts found for " + inputSam + ".  Your reads may not be paired properly.");
        }
	}
	
	/**
	 * Convert the input SAM/BAM file into a single fastq file.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 */
	public void convert(String inputSam, String outputFastq) throws IOException {
		String last1Read = "";
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);

        output1 = new FastqOutputFile();
        output1.init(outputFastq);
        int lineCnt = 0;
        
        for (SAMRecord read : reader) {
    		if (!read.getReadName().equals(last1Read)) {
    			output1.write(samReadToFastqRecord(read));
    			last1Read = read.getReadName();
    		}
    		
            lineCnt++;
            if ((lineCnt % 1000000) == 0) {
                System.out.println("record: " + lineCnt);
            }
        }
                
        output1.close();
        reader.close();
	}
	
	private FastqRecord samReadToFastqRecord(SAMRecord read) {
		String bases = read.getReadString();
		String qualities = read.getBaseQualityString();
		
		if ((isMapspliceFusions) && (isFusion(read))) {
			String fusionAttribute = (String) read.getAttribute("ZF");
			int parenIdx = fusionAttribute.indexOf('(');
			if ((parenIdx < 0) || (parenIdx+2 >= fusionAttribute.length())) {
				throw new IllegalArgumentException("Invalid fusion attribute for read: [" + read.getSAMString() + "]");
			}
			boolean isDonerNegative = fusionAttribute.charAt(parenIdx+1) == '-';
			boolean isAccepterNegative = fusionAttribute.charAt(parenIdx+2) == '-';
			
			int donerLength = 0;
			int accepterLength = 0;
			
			CigarElement firstElement = read.getCigar().getCigarElement(0);
			CigarElement lastElement  = read.getCigar().getCigarElement(read.getCigar().getCigarElements().size()-1);
			
			if (firstElement.getOperator() == CigarOperator.S) {
				donerLength = firstElement.getLength();
				accepterLength = read.getReadLength() - donerLength;
			} else if (lastElement.getOperator() == CigarOperator.S) {
				accepterLength = lastElement.getLength();
				donerLength = read.getReadLength() - accepterLength;
			} else {
				throw new IllegalArgumentException (
						"Expected soft clip at begin or end of Cigar for fusion read: [" + read.getSAMString() + "]");
			}
			
			String donerBases = read.getReadString().substring(0, donerLength);
			String accepterBases = read.getReadString().substring(donerLength, read.getReadLength());
			
			if (isDonerNegative) {
				donerBases = reverseComplementor.reverseComplement(donerBases);
			}
			
			if (isAccepterNegative) {
				accepterBases = reverseComplementor.reverseComplement(accepterBases);
			}
			
			bases = donerBases + accepterBases;
		}
		else if (read.getReadNegativeStrandFlag()) {
			bases = reverseComplementor.reverseComplement(bases);
			qualities = reverseComplementor.reverse(qualities);
		}
		
		FastqRecord fastq = new FastqRecord("@" + read.getReadName(), bases, qualities);
		
		return fastq;
	}
	
	private boolean isFusion(SAMRecord read) {
		return read.getAttribute("ZF") != null;
	}
	
	private boolean isFirstInPair(SAMRecord read) {
		boolean isFirstInPair;
		
		if (shouldIdentifyEndByReadId) {
			isFirstInPair = read.getReadName().endsWith(end1Suffix);
			
		} else {
			isFirstInPair = read.getFirstOfPairFlag();
		}
		
		return isFirstInPair;
	}
	
	private boolean isSecondInPair(SAMRecord read) {
		boolean isSecondInPair;
		
		if (shouldIdentifyEndByReadId) {
			isSecondInPair = read.getReadName().endsWith(end2Suffix);
			
		} else {
			isSecondInPair = read.getSecondOfPairFlag();
		}
		
		return isSecondInPair;
	}
	
	public void setEndSuffixes(String end1Suffix, String end2Suffix) {
		this.shouldIdentifyEndByReadId = true;
		this.end1Suffix = end1Suffix;
		this.end2Suffix = end2Suffix;
	}
	
	public void setMapspliceFusions(boolean isMapspliceFusions) {
		this.isMapspliceFusions = isMapspliceFusions;
	}

/*
	public static void run(String[] args) throws IOException {
		Sam2FastqOptions options = new Sam2FastqOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			long s = System.currentTimeMillis();
			System.out.println("sam2fastq starting");
			
			Sam2Fastq sam2Fastq = new Sam2Fastq();
			if (options.isPairedEnd()) {
				
				if (options.shouldIdEndByReadName()) {
					sam2Fastq.setEndSuffixes(options.getEnd1Suffix(), options.getEnd2Suffix());
				}
				
				sam2Fastq.setMapspliceFusions(options.isMapspliceFusions());
				
				sam2Fastq.convert(options.getInputFile(), options.getFastq1(), options.getFastq2());
			} else {
				sam2Fastq.convert(options.getInputFile(), options.getFastq1());
			}
			
			long e = System.currentTimeMillis();
			System.out.println("sam2fastq done.  Elapsed secs: " + (e-s)/1000);
		}
	}
*/
/*	
	public static void main(String[] args) throws Exception {
		String argz = "--end1 /1 --end2 /2 --in /home/lisle/sam2fastq/round2/test2.sam  --fastq1 /home/lisle/sam2fastq/round2/1.fastq --fastq2 /home/lisle/sam2fastq/round2/2.fastq --mapsplice";
		
		run(argz.split(" "));
	}
*/
}
