/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 * Utility class for assessing mismatches by base quality
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class MismatchBaseQuality {
	
	//TODO: Move to utils
	private CompareToReference2 c2r;
	private Map<Integer, Long> qualityCounts;
	private Map<Integer, Long> qualityMismatches;

	public void run(String input, String ref) throws IOException {
		c2r = new CompareToReference2();
		c2r.init(ref);
		
		qualityMismatches = new HashMap<Integer, Long>();
		qualityCounts = new HashMap<Integer, Long>();
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		
		int min = 129;
		int max = -129;
		
		for (SAMRecord read : reader) {
			if (!read.getReadUnmappedFlag()) {
				List<Integer> mismatchPositions = c2r.mismatchPositions(read);
				
				for (int pos : mismatchPositions) {
					int qual = read.getBaseQualityString().charAt(pos) - 33;
					if (qual > max) {
						max = qual;
					}
					
					if (qual < min) {
						min = qual;
					}
					
					addMismatch(qual);
				}
				
				for (int i=0; i< read.getBaseQualityString().length(); i++) {
					int qual = read.getBaseQualityString().charAt(i) - 33;
					
					if (qual > max) {
						max = qual;
					}
					
					if (qual < min) {
						min = qual;
					}
					
					addCount(qual);
				}
			}
		}
		
		reader.close();
		
		for (int i=min; i<max; i++) {
			long count = getValue(qualityCounts, i);
			long mismatches = getValue(qualityMismatches, i);
			System.out.println(i + "," + count + "," +
					mismatches + "," + ((double) mismatches / (double) count));
		}
	}
	
	private long getValue(Map<Integer, Long> map, int pos) {
		Long l = map.get(pos);
		if (l != null) { 
			return l.longValue();
		} else {
			return 0L;
		}
	}

	private void addCount(int pos) {
		Long count = qualityCounts.get(pos);
		if (count == null) {
			qualityCounts.put(pos, 1L);
		} else {
			qualityCounts.put(pos, count + 1);
		}
	}
	
	private void addMismatch(int pos) {
		Long count = qualityMismatches.get(pos);
		if (count == null) {
			qualityMismatches.put(pos, 1L);
		} else {
			qualityMismatches.put(pos, count + 1);
		}
	}
	
	public static void main(String[] args) throws Exception {
		new MismatchBaseQuality().run(args[0], args[1]);
	}
}
