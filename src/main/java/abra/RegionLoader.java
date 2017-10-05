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

	public List<Feature> load(String regionFile, boolean hasPresetKmers) throws FileNotFoundException, IOException {
		List<Feature> features = new ArrayList<Feature>();
		
		int start = BED_START_IDX;
		int end = BED_END_IDX;
				
		BufferedReader reader = new BufferedReader(new FileReader(regionFile));
		
		try {
			
			String lastChr = "";
			long lastStart = -1;
			
			String line = reader.readLine();
			
			int cnt = 0;
			
			while (line != null) {
				if(line.startsWith("#") || line.trim().isEmpty() || line.startsWith("track") || line.startsWith("browser")) {
					line = reader.readLine();
					continue;
				}
				String[] fields = line.split("\t");
				
				String chromosome = fields[SEQNAME_IDX];
				long startPos = Long.valueOf(fields[start]);
				long endPos = Long.valueOf(fields[end]);
				
				if (startPos > endPos) {
					throw new IllegalArgumentException("Region end must be greater than region start in target BED file: " + line);
				}
				
				if (lastChr.equals(chromosome) && startPos < lastStart) {
					throw new IllegalArgumentException("Target BED file must be sorted in increasing coordinate order (grouped by chromosome): " + line);
				}
				
				Feature feature = new Feature(chromosome, startPos, endPos);
				
				if (fields.length >= KMER_SIZE_IDX+1 && hasPresetKmers) {
					int kmerSize = Integer.parseInt(fields[KMER_SIZE_IDX]);
					feature.setKmer(kmerSize);
				}
				
				features.add(feature);
				
				line = reader.readLine();
				cnt++;
				
				lastChr = chromosome;
				lastStart = startPos;
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
		
		Logger.info("Collapsed regions from " + regions.size() + " to " + collapsedRegions.size());
		
		return collapsedRegions;
	}
	
	public static void main(String[] args) throws Exception {
		RegionLoader loader = new RegionLoader();
//		List<Feature> regions = loader.load("/home/lmose/dev/abra/issue12/test.bed");
		List<Feature> regions = loader.load("/home/lmose/dev/abra/issue12/test2.bed", false);
		
		regions = RegionLoader.collapseRegions(regions, 100);

		/*
		for (Feature region : regions) {
			if (region.getLength() <= 0) {
				System.out.println(region + " - " + region.getLength());	
			}
			
		}
		*/

		regions = ReAligner.splitRegions(regions);	
		
		for (Feature region : regions) {
			if (region.getLength() <= 0) {
				System.err.println(region + " - " + region.getLength());	
			}
			
		}
	}
}
