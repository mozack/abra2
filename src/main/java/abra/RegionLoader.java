/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads region info into memory.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class RegionLoader {
	
	private static final int SEQNAME_IDX = 0;
	private static final int GTF_START_IDX   = 3;
	private static final int GTF_END_IDX     = 4;
	private static final int BED_START_IDX = 1;
	private static final int BED_END_IDX = 2;

	public List<Feature> load(String regionFile) throws FileNotFoundException, IOException {
		List<Feature> features = new ArrayList<Feature>();
		
		// Assume BED format unless extension is .gtf
		int start = BED_START_IDX;
		int end = BED_END_IDX;
		
		if (regionFile.endsWith(".gtf")) {
			start = GTF_START_IDX;
			end = GTF_END_IDX;
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(regionFile));
		
		String line = reader.readLine();
		
		while (line != null) {
			String[] fields = line.split("\t");			
			features.add(new Feature(fields[SEQNAME_IDX], Long.valueOf(fields[start]), Long.valueOf(fields[end])));
			
			line = reader.readLine();
		}
		
		reader.close();
		
		return features;
	}
}
