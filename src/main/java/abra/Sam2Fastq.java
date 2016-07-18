/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.io.IOException;
import java.util.List;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

/**
 * Responsible for creating fastq files for reads of interest to be realigned.
 * Those reads that are not eligible for realignment are written directly to the output BAM file.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Sam2Fastq {
	
	public static final String FIELD_DELIMITER = "~|";
	private static final int MAX_SAM_READ_NAME_LENGTH = 255;
	
	public static final int MIN_OFF_TARGET_MAPQ = 30;
	
	public static int COMPRESSION_LEVEL = 1;
	
	private FastqOutputFile output1;
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	private boolean shouldIdentifyEndByReadId = false;
	private String end1Suffix;
	private String end2Suffix;
	private RegionTracker regionTracker;
		
	/**
	 * Convert the input SAM/BAM file into a single fastq file.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 */
	public void convert(String inputSam, String outputFile, CompareToReference2 c2r,
			SAMFileWriter writer, boolean isPairedEnd,
			List<Feature> regions, int minMappingQuality, boolean isBamOutput) throws IOException {
		
		System.err.println("sam: " + inputSam);
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);
        
        SAMFileWriter toProcessWriter = null;
        
        if (isBamOutput) {
        
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
			SAMFileHeader header = reader.getFileHeader();
			header.setSortOrder(SortOrder.unsorted);
			
			toProcessWriter = writerFactory.makeBAMWriter(
					header, false, new File(outputFile), COMPRESSION_LEVEL);
		
        } else {
	        output1 = new FastqOutputFile();
	        output1.init(outputFile);
        }
        
        int lineCnt = 0;
        
        this.regionTracker = new RegionTracker(regions, reader.getFileHeader());
        
        for (SAMRecord read : reader) {
        	if (!SAMRecordUtils.isPrimary(read)) {
        		// Write secondary / supplemental reads directly to the output BAM file.
        		// TODO: If primary is realigned, perhaps this should be handled differently?
        		
        		if (writer != null) {
        			writer.addAlignment(read);
        		}
        	}
        	else if (!SAMRecordUtils.isFiltered(isPairedEnd, read)) {
    			
    			// These tags can be lengthy, so remove them.
    			// TODO: Improve the way this is handled
    			read.setAttribute("XA", null);
    			read.setAttribute("OQ", null);
    			read.setAttribute("MD", null);
    			read.setAttribute("BQ", null);
    			read.setAttribute("BI", null);
    			read.setAttribute("BD", null);
    			
    			int yx = 0;
    			
    			boolean isAmbiguous = !read.getReadUnmappedFlag() && read.getMappingQuality() == 0;
    			
    			if ((!read.getReadFailsVendorQualityCheckFlag()) && (!isAmbiguous)) {
	    			// Calculate the number of mismatches to reference for this read.
	    			if (c2r != null) {
	    				try {
	    					yx = SAMRecordUtils.getEditDistance(read, c2r);
	    				} catch (ArrayIndexOutOfBoundsException e) {
	    					System.err.println("Index error for read: " + read.getSAMString());
	    					throw e;
	    				}
	    			} else {
	    				yx = read.getReadLength();
	    			}
	    			
	    			read.setAttribute("YX", yx);
    			}

    			boolean offTargetFiltered = false;
				if (yx > 0 && !read.getReadUnmappedFlag() && read.getMappingQuality() < MIN_OFF_TARGET_MAPQ && !regionTracker.isInRegion(read)) {
					read.setAttribute("YR", 2);
					offTargetFiltered = true;
//					System.out.println("Filtering off target: " + read.getSAMString());
				}
    			
    			if ((yx > 0 && !offTargetFiltered && read.getMappingQuality() >= minMappingQuality) || (read.getReadUnmappedFlag())) {
    				
	    			if ((!read.getReadUnmappedFlag()) && (!regionTracker.isInRegion(read))) {
	    				read.setAttribute("YR", 1);
	    			}
    					    			
	    			// SA tag causes read info to not fit into read name, remove it for now.
	    			read.setAttribute("SA", null);
	    			
	    			try {
	    				if (isBamOutput) {
	    					toProcessWriter.addAlignment(samReadToUnmappedSam(read));
	    				} else {
	    					output1.write(samReadToFastqRecord(read));
	    				}
	    			} catch (IllegalArgumentException e) {
	    				System.err.println("Error on: " + read.getSAMString());
	    				e.printStackTrace();
	    				throw e;
	    			}
    			} else if (writer != null) {
    				// Either xactly matches reference or failed vendor QC or low mapq, so
    				// output directly to final BAM
    				writer.addAlignment(read);
    			}
    		}
    		
            lineCnt++;
            if ((lineCnt % 1000000) == 0) {
                System.err.println("record: " + lineCnt);
            }
        }
                
        if (isBamOutput) {
        	toProcessWriter.close();
        } else {
        	output1.close();
        }
        reader.close();
	}

	private SAMRecord samReadToUnmappedSam(SAMRecord read) {
		
		String bases = read.getReadString();
		String qualities = read.getBaseQualityString();
		
		if (read.getReadNegativeStrandFlag()) {
			bases = reverseComplementor.reverseComplement(bases);
			qualities = reverseComplementor.reverse(qualities);
		}
		
		read.setReadString("");
		read.setBaseQualityString("");
		
		String readStr = read.getSAMString();

		readStr = readStr.replace("\t", FIELD_DELIMITER);
		
		// Trim newline if applicable
		if (readStr.charAt(readStr.length()-1) == '\n') {
			readStr = readStr.substring(0, readStr.length()-1);
		}
		
		String readName = readStr;
		
		read.setReadName(readName);
		read.setReadString(bases);
		read.setBaseQualityString(qualities);
		read.setFlags(0);
		read.clearAttributes();
		read.setReadUnmappedFlag(true);
		read.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);
		read.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
		read.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
		read.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);

		return read;
	}
	
	private FastqRecord samReadToFastqRecord(SAMRecord read) {
		
		String bases = read.getReadString();
		String qualities = read.getBaseQualityString();
		
		if (read.getReadNegativeStrandFlag()) {
			bases = reverseComplementor.reverseComplement(bases);
			qualities = reverseComplementor.reverse(qualities);
		}
		
		read.setReadString("");
		read.setBaseQualityString("");
		
		String readStr = read.getSAMString();

		readStr = readStr.replace("\t", FIELD_DELIMITER);
		
		// Trim newline if applicable
		if (readStr.charAt(readStr.length()-1) == '\n') {
			readStr = readStr.substring(0, readStr.length()-1);
		}
		
		String readName = readStr; 
		
		FastqRecord fastq = new FastqRecord("@" + readName, bases, qualities);
		
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

	public static void main(String[] args) throws Exception {

		/*
		String inputSam = "/home/lmose/dev/abra/s2fq_test/chr22.bam";
		String outputSam = "/home/lmose/dev/abra/s2fq_test/chr22_out.bam";
		String bed = "/home/lmose/dev/abra/s2fq_test/chr22.bed";
		//String tempReadFile = "/home/lmose/dev/abra/s2fq_test/chr22.fastq.gz";
		String tempReadFile = "/home/lmose/dev/abra/s2fq_test/chr22.temp.bam";
		String ref = "/home/lmose/reference/chr22/chr22.fa";
		*/
		
		String inputSam = args[0];
		String outputSam = args[1];
		String bed = args[2];
		String tempReadFile = args[3];
		String ref = args[4];
		COMPRESSION_LEVEL  = Integer.parseInt(args[5]);
		
		
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		SAMFileHeader header = reader.getFileHeader();
		reader.close();
				
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		header.setSortOrder(SortOrder.unsorted);
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				header, false, new File(outputSam));
		
		Sam2Fastq s2f = new Sam2Fastq();
		
		
		RegionLoader loader = new RegionLoader();
		List<Feature> regions = loader.load(bed, false);
				
		long s = System.currentTimeMillis();
		s2f.convert(inputSam, tempReadFile, c2r, writer, false, 
				regions, 20, true);
		long e = System.currentTimeMillis();
		
		
		writer.close();
		
		System.err.println("Elapsed: " + (e-s)/1000);
	}
}
