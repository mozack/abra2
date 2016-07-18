package abra;

import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;

public class SVReadCounter {
	
	private Map<String, Integer> breakpointCounts = new HashMap<String, Integer>();
	
	private static final int MAX_EDIT_DISTANCE = 5;
	
	private SAMLineParser parser;
	
	private Map<String, Integer> counts;

	public Map<String, Integer> countReadsSupportingBreakpoints(SamReader reader, int readLength, SAMFileHeader samHeader) {
		
		parser = new SAMLineParser(new DefaultSAMRecordFactory(),
                ValidationStringency.SILENT, samHeader,
                null, null);
		
		String fullMatch = readLength + "M";
		
		// Require 90% of the read to overlap the breakpoint
		int minStart = (int) (readLength * .10);
		int maxStart = (int) (readLength *.9) + 1;
		
		// TODO: Need way to query mapped reads only
		for (SAMRecord read : reader) {
			if (!read.getReadUnmappedFlag() && read.getCigarString().equals(fullMatch)) {
				if (read.getAlignmentStart() >= minStart && read.getAlignmentStart() <= maxStart) {
					int editDistance = SAMRecordUtils.getIntAttribute(read, "NM");
					
					if (editDistance <= MAX_EDIT_DISTANCE) {
						SAMRecord orig = getOrigRecord(read, samHeader);
						int origEditDistance = SAMRecordUtils.getIntAttribute(orig, "YX");
						if (editDistance < origEditDistance) {
							//TODO: Inspect alternate alignments
							/*
							String[] regions = read.getReferenceName().split("__");
							if (regions.length == 2) {
								String[] region1 = regions[0].split("_");
								String[] region2 = regions[1].split("_");
								if (region1.length >= 3 && region2.length >= 5) {
									int region1IdxPad = region1.length - 3; // Adjustment for underscore in chromosome name
									int region2IdxPad = region2.length - 5; // Adjustment for underscore in chromosome name
									
									String chr1 = region1[0];
									for (int i=1; i<=region1IdxPad; i++) {
										chr1 += "_";
										chr1 += region1[i];
									}
									int start1 = Integer.parseInt(region1[1 + region1IdxPad]);
									int stop1 = Integer.parseInt(region1[2 + region1IdxPad]);
									
									String chr2 = region2[0];
									for (int i=1; i<=region2IdxPad; i++) {
										chr1 += "_";
										chr2 += region2[i];
									}
									int start2 = Integer.parseInt(region2[1 + region2IdxPad]);
									int stop2 = Integer.parseInt(region2[2 + region2IdxPad]);
								}
							}
							*/
							
							String[] refFields = read.getReferenceName().split("_");
							if (refFields.length >= 8) {
//								String breakpointGroupId = refFields[0] + "_" + refFields[1] + "\t" + refFields[2] + ":" + refFields[3] + "\t" +
//										refFields[4] + ":" + refFields[5];

								
								String breakpointGroupId = refFields[0] + "_" + refFields[1] + "\t" + refFields[2] + ":" + refFields[3] + "\t" +
										refFields[4] + ":" + refFields[5] + "\t" + refFields[refFields.length-2] + "\t" + refFields[refFields.length-1];
								

								Integer count = breakpointCounts.get(breakpointGroupId);
								if (count == null) {
									breakpointCounts.put(breakpointGroupId, 1);
								} else {
									breakpointCounts.put(breakpointGroupId, count + 1);
								}
							} else {
								System.err.println("Error analyzing breakpoint for: " + read.getSAMString());
							}
						}
					}
				}
			}
		}
		
		this.counts = breakpointCounts;
		
		return breakpointCounts;
	}
	
	private SAMRecord getOrigRecord(SAMRecord read, SAMFileHeader samHeader) {
		String origSamStr = read.getReadName();
		origSamStr = origSamStr.replace(Sam2Fastq.FIELD_DELIMITER, "\t");
		SAMRecord orig;
		try {
			orig = parser.parseLine(origSamStr);
		} catch (RuntimeException e) {
			System.err.println("Error processing: [" + origSamStr + "]");
			System.err.println("Contig read: [" + read.getSAMString() + "]");
			e.printStackTrace();
			throw e;
		}
		orig.setHeader(samHeader);

		return orig;
	}
	
	public Map<String, Integer> getCounts() {
		return counts;
	}
	
	public static void main(String[] args) throws Exception {
		String file = "/home/lmose/dev/abra/sv/virus_test2/t.sv.bam";
		
		final SamReader reader =
		        SamReaderFactory.make()
		                .validationStringency(ValidationStringency.SILENT)
		                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
		                .open(SamInputResource.of(file));
		
		SAMFileHeader header = reader.getFileHeader();
		
		System.err.println("header: " + header);

		SVReadCounter counter = new SVReadCounter();
		Map<String, Integer> counts = counter.countReadsSupportingBreakpoints(reader, 100, header);
		reader.close();
		
		System.err.println(counts);
	}
}
