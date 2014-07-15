/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.io.IOException;
import java.util.List;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Responsible for creating fastq files for reads of interest to be realigned.
 * Those reads that are not eligible for realignment are written directly to the output BAM file.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class Sam2Fastq {
	
	public static final String FIELD_DELIMITER = "~|";
	private static final int MAX_SAM_READ_NAME_LENGTH = 255;
	
	private FastqOutputFile output1;
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	private boolean shouldIdentifyEndByReadId = false;
	private boolean isMapspliceFusions = false;
	private String end1Suffix;
	private String end2Suffix;
	private RegionTracker regionTracker;
		
	/**
	 * Convert the input SAM/BAM file into a single fastq file.
	 * Input SAM files that contain multiple mappings should be sorted by read name.
	 */
	public void convert(String inputSam, String outputFastq, CompareToReference2 c2r,
			SAMFileWriter writer, boolean isPairedEnd,
			List<Feature> regions) throws IOException {
		
		System.out.println("sam: " + inputSam);
		
        SAMFileReader reader = new SAMFileReader(new File(inputSam));
        reader.setValidationStringency(ValidationStringency.SILENT);
        
        output1 = new FastqOutputFile();
        output1.init(outputFastq);
        int lineCnt = 0;
        
        this.regionTracker = new RegionTracker(regions, reader.getFileHeader());
        
        for (SAMRecord read : reader) {
    		if (SAMRecordUtils.isPrimary(read) && (!SAMRecordUtils.isFiltered(isPairedEnd, read))) {
    			
    			// These tags can be lengthy, so remove them.
    			// TODO: Improve the way this is handled
    			read.setAttribute("XA", null);
    			read.setAttribute("OQ", null);
    			read.setAttribute("MD", null);
    			read.setAttribute("BQ", null);
    			
    			int yx = 0;
    			
//    			boolean isAmbiguous = !read.getReadUnmappedFlag() && read.getMappingQuality() == 0;
    			
    			boolean isAmbiguous = false;
    			
    			if ((!read.getReadFailsVendorQualityCheckFlag()) && (!isAmbiguous)) {
	    			// Calculate the number of mismatches to reference for this read.
	    			if (c2r != null) {
	    				yx = SAMRecordUtils.getEditDistance(read, c2r);
	    			} else {
	    				yx = read.getReadLength();
	    			}
	    			
	    			read.setAttribute("YX", yx);
    			}
    			
    			if (yx > 0) {
    				
//	    			if ((!read.getReadUnmappedFlag()) && (!regionTracker.isInRegion(read))) {
//	    				read.setAttribute("YR", 1);
//	    			}
    				
    				// Does not exactly match reference, output FASTQ record
	    			output1.write(samReadToFastqRecord(read, c2r));
    			} else if (writer != null) {
    				// Either xactly matches reference or failed vendor QC, so
    				// output directly to final BAM
    				writer.addAlignment(read);
    			}
    		}
    		
            lineCnt++;
            if ((lineCnt % 1000000) == 0) {
                System.out.println("record: " + lineCnt);
            }
        }
                
        output1.close();
        reader.close();
	}
	
	private FastqRecord samReadToFastqRecord(SAMRecord read, CompareToReference2 c2r) {
		
		String bases = read.getReadString();
		String qualities = read.getBaseQualityString();
		
		if (read.getReadNegativeStrandFlag()) {
			bases = reverseComplementor.reverseComplement(bases);
			qualities = reverseComplementor.reverse(qualities);
		}
		
		read.setReadString("");
		read.setBaseQualityString("");
		
		String readStr = read.getSAMString();
/*		
		if (readStr.length() > MAX_SAM_READ_NAME_LENGTH) {
			String msg = "Warning!  Max SAM Read name length exceeded for: " + readStr;
			System.out.println(msg);
//			throw new RuntimeException(msg);
		}
*/

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

	public static void main(String[] args) throws Exception {
		
		String inputSam = "/home/lmose/dev/ayc/opt/opt2/s2f/t1.bam";
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		SAMFileHeader header = reader.getFileHeader();
		reader.close();
				
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("/home/lmose/reference/chr1/chr1.fa");
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
//		writerFactory.setUseAsyncIo(true);
		
		header.setSortOrder(SortOrder.unsorted);
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				header, false, new File("/home/lmose/dev/ayc/opt/opt2/s2f/output_t1.bam"));
		
		Sam2Fastq s2f = new Sam2Fastq();
		
		
		RegionLoader loader = new RegionLoader();
		List<Feature> regions = loader.load("/home/lmose/dev/ayc/regions/clinseq5/uncseq5.gtf");
		
		regions = RegionLoader.collapseRegions(regions, 100);
		
		regions = ReAligner.splitRegions(regions);		
		
		long s = System.currentTimeMillis();
		s2f.convert(inputSam, "/home/lmose/dev/ayc/opt/opt2/s2f/t2.fastq.gz", c2r, writer, false, 
				regions);
		long e = System.currentTimeMillis();
		
		
		writer.close();
		
		System.out.println("Elapsed: " + (e-s)/1000);
		
		
//		s2f.convert("/home/lmose/dev/abra_wxs/21_1071/small_tumor.abra.bam", "/home/lmose/dev/abra_wxs/21_1071/t.fastq", c2r);
	}
	/*
	public static void main(String[] args) throws Exception {
		
		String inputSam = "/home/lmose/dev/ayc/opt/t7.bam";
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		SAMFileHeader header = reader.getFileHeader();
		reader.close();
				
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("/home/lmose/reference/chr7/chr7.fa");
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
//		writerFactory.setUseAsyncIo(true);
		
		header.setSortOrder(SortOrder.unsorted);
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				header, false, new File("/home/lmose/dev/ayc/opt/output_t7.bam"));
		
		Sam2Fastq s2f = new Sam2Fastq();
		
		long s = System.currentTimeMillis();
		
		
		System.out.println("chr1: " + header.getSequenceIndex("chr1"));
		System.out.println("chr2: " + header.getSequenceIndex("chr2"));
		System.out.println("chrY: " + header.getSequenceIndex("chrY"));
		System.out.println("chrM: " + header.getSequenceIndex("chrM"));
		
//		String inputSam, String outputFastq, CompareToReference2 c2r,
//		SAMFileHeader header, SAMFileWriter writer, ReAligner realigner
		
		GtfLoader loader = new GtfLoader();
		List<Feature> regions = loader.load("/home/lmose/dev/ayc/regions/clinseq5/uncseq5.gtf");
		
		regions = ReAligner.collapseRegions(regions, 100);
		
		regions = ReAligner.splitRegions(regions);		
		
		s2f.convert(inputSam, "/home/lmose/dev/ayc/opt/t7.fastq.gz", c2r, header, writer, false, 
				regions);
		long e = System.currentTimeMillis();
		
		System.out.println("Elapsed: " + (e-s)/1000);
		
		writer.close();
		
//		s2f.convert("/home/lmose/dev/abra_wxs/21_1071/small_tumor.abra.bam", "/home/lmose/dev/abra_wxs/21_1071/t.fastq", c2r);
	}
	*/
}
