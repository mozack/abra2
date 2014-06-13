/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * Loads region info into memory.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class RegionLoader {
	
	private static final int SEQNAME_IDX = 0;
	private static final int BED_START_IDX = 1;
	private static final int BED_END_IDX = 2;
	private static final int KMER_SIZE_IDX = 3;

	public List<Feature> load(String regionFile, boolean requireKmerSizes) throws FileNotFoundException, IOException {
		List<Feature> features = new ArrayList<Feature>();
		
		int start = BED_START_IDX;
		int end = BED_END_IDX;
				
		BufferedReader reader = new BufferedReader(new FileReader(regionFile));
		
		try {
			String line = reader.readLine();
			
			int cnt = 0;
			
			while (line != null) {
				String[] fields = line.split("\t");
				Feature feature = new Feature(fields[SEQNAME_IDX], Long.valueOf(fields[start]), Long.valueOf(fields[end])); 
				
				if (fields.length >= KMER_SIZE_IDX+1) {
					int kmerSize = Integer.parseInt(fields[KMER_SIZE_IDX]);
					feature.setKmer(kmerSize);
				} else if (requireKmerSizes) {
					throw new IllegalArgumentException("Kmer size missing in bed file.  Region: [" + feature.getDescriptor() + "].  " +
							"Please use KmerSizeEvaluator to generate kmer sizes for regions or specify kmer sizes using the --kmer option");
				}
				
				features.add(feature);
				
				line = reader.readLine();
				cnt++;
				if ((cnt % 100000) == 0) {
					System.out.println("Loaded " + cnt + " regions");
					System.out.flush();
				}
			}
		} finally {
			reader.close();
		}
		
		return features;
	}
	
	public static List<Feature> collapseRegions(Collection<Feature> regions, int maxGap) {
		List<Feature> collapsedRegions = new ArrayList<Feature>();
		
		Feature currentRegion = null;
		
		for (Feature region : regions) {
			if (currentRegion != null) {
				if ((currentRegion.getSeqname().equals(region.getSeqname())) && 
					(currentRegion.getEnd() + (maxGap) >= region.getStart())) {
					
					currentRegion.setEnd(region.getEnd());
				} else {
					collapsedRegions.add(currentRegion);
					currentRegion = region;
				}
			} else {
				currentRegion = region;
			}
		}
		
		if (currentRegion != null) {
			collapsedRegions.add(currentRegion);
		}
		
		System.out.println("Collapsed regions from " + regions.size() + " to " + collapsedRegions.size());
		
		return collapsedRegions;
	}
}
