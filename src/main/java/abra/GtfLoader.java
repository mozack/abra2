/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Loads GTF file content into memory.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class GtfLoader {
	
	private static final int SEQNAME_IDX = 0;
	private static final int START_IDX   = 3;
	private static final int END_IDX     = 4;

	public List<Feature> load(String gtfFile) throws FileNotFoundException, IOException {
		List<Feature> features = new ArrayList<Feature>();
		
		BufferedReader reader = new BufferedReader(new FileReader(gtfFile));
		
		String line = reader.readLine();
		
		while (line != null) {
			String[] fields = line.split("\t");			
			features.add(new Feature(fields[SEQNAME_IDX], Long.valueOf(fields[START_IDX]), Long.valueOf(fields[END_IDX])));
			
			line = reader.readLine();
		}
		
		reader.close();
		
		return features;
	}	
}
