package abra;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;

public class SVReadCounter {
	
	private Map<String, Integer> breakpointCounts = new HashMap<String, Integer>();
	
	private static final int MAX_EDIT_DISTANCE = 5;
	
	private SAMLineParser parser;

	public Map<String, Integer> countReadsSupportingBreakpoints(String svSam, int readLength, SAMFileHeader samHeader) {
		
		parser = new SAMLineParser(samHeader);
		
		SAMFileReader reader = new SAMFileReader(new File(svSam));
		reader.setValidationStringency(ValidationStringency.SILENT);

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
							String[] refFields = read.getReferenceName().split("_");
							if (refFields.length >= 6) {
								String breakpointGroupId = refFields[0] + "_" + refFields[1] + "\t" + refFields[2] + ":" + refFields[3] + "\t" +
										refFields[4] + ":" + refFields[5];
								Integer count = breakpointCounts.get(breakpointGroupId);
								if (count == null) {
									breakpointCounts.put(breakpointGroupId, 1);
								} else {
									breakpointCounts.put(breakpointGroupId, count + 1);
								}
							} else {
								System.out.println("Error analyzing breakpoint for: " + read.getSAMString());
							}
						}
					}
				}
			}
		}
		
		reader.close();
		
		return breakpointCounts;
	}
	
	private SAMRecord getOrigRecord(SAMRecord read, SAMFileHeader samHeader) {
		String origSamStr = read.getReadName();
		origSamStr = origSamStr.replace(Sam2Fastq.FIELD_DELIMITER, "\t");
		SAMRecord orig;
		try {
			orig = parser.parseLine(origSamStr);
		} catch (RuntimeException e) {
			System.out.println("Error processing: [" + origSamStr + "]");
			System.out.println("Contig read: [" + read.getSAMString() + "]");
			e.printStackTrace();
			throw e;
		}
		orig.setHeader(samHeader);
		
//		orig.setReadString(read.getReadString());
//		orig.setBaseQualityString(read.getBaseQualityString());

		return orig;
	}
	
	public static void main(String[] args) {
		SVReadCounter svc = new SVReadCounter();
		
		SAMFileReader reader = new SAMFileReader(new File("/home/lmose/dev/abra/sv/test_sv.sam"));
		
		Map<String, Integer> counts = svc.countReadsSupportingBreakpoints("/home/lmose/dev/abra/sv/sv.bam", 100, reader.getFileHeader());
		
		for (String bp : counts.keySet()) {
			System.out.println(bp + "\t" + counts.get(bp));
		}
		
		reader.close();
	}
}
