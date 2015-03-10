package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

/**
 * Keeps track of current region and provides
 * utility methods on reads related to regions.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class RegionTracker {
	
	private List<Feature> regions;
	private Feature currentRegion;
	private Iterator<Feature> regionIter;
	private SAMFileHeader header;

	public RegionTracker(List<Feature> regions, SAMFileHeader header) {
		this.regions = regions;
		this.header = header;
		init();
	}
	
	public void init() {
		regionIter = regions.iterator();
        if (regionIter.hasNext()) {
        	currentRegion = regionIter.next();
        }		
	}
	
	public boolean isInRegion(SAMRecord read) {
		
		while ((currentRegion != null) &&
			   (!currentRegion.overlapsRead(read)) &&
			   (!isRegionBeyondRead(header, currentRegion, read))) {
			
			if (regionIter.hasNext()) {
				currentRegion = (Feature) regionIter.next();
			} else {
				currentRegion = null;
			}
		}
		
		return (currentRegion != null) && (currentRegion.overlapsRead(read));
	}
	
	public List<Feature> identifyTargetRegions(List<String> files, int minBaseQuality, 
			int readLength, CompareToReference2 c2r) {
		
		Map<String, List<Integer>> locations = new HashMap<String, List<Integer>>();
		
		for (String file : files) {
	        SAMFileReader reader = new SAMFileReader(new File(file));
	        reader.setValidationStringency(ValidationStringency.SILENT);

	        header = reader.getFileHeader();
	        init();
	        
	        for (SAMRecord read : reader) {
	        	if (isInRegion(read)) {
	        		if (read.getCigarString().contains("I") || read.getCigarString().contains("D")) {
	        			addLocation(read, locations);
	        		} else if (read.getCigarString().contains("S") && c2r.numHighQualityMismatches(read, minBaseQuality) > 1) {
	        			addLocation(read, locations);
	        		}
	        	}
	        }
	        
	        reader.close();
		}
		
		return locationsToRegions(locations, readLength);
	}
	
	private List<Feature> locationsToRegions(Map<String, List<Integer>> locations, int readLength) {
		List<Feature> regions = new ArrayList<Feature>();
		
		for (String chr : locations.keySet()) {
			int prev = -readLength;
			int start = -readLength;
			
			List<Integer> positions = locations.get(chr);
			
			for (int pos : positions) {
				
				if (start < 0) {
					start = pos;
				}
				
				if (pos < prev+readLength) {
					prev = pos;
				} else {
					if (start > 0) {
						int end = prev + 2*readLength;
						start = start-readLength;
						regions.add(new Feature(chr, start, end));
						start = pos;
						prev = pos;
					}
				}
			}
		}
		
		regions = RegionLoader.collapseRegions(regions, readLength);
		regions = ReAligner.splitRegions(regions);
		
		return regions;
	}
	
	private void addLocation(SAMRecord read, Map<String, List<Integer>> locations) {
		if (!locations.containsKey(read.getReferenceName())) {
			locations.put(read.getReferenceName(), new ArrayList<Integer>());
		}
		
		locations.get(read.getReferenceName()).add(read.getAlignmentStart());
	}
	
	private boolean isRegionBeyondRead(SAMFileHeader header, Feature region, SAMRecord read) {
		boolean isRegionBeyond = false;
		
		int regionChrIdx = header.getSequenceIndex(region.getSeqname());
		int readChrIdx = header.getSequenceIndex(read.getReferenceName());
		
		if (regionChrIdx > readChrIdx) {
			isRegionBeyond = true;
		} else if (regionChrIdx < readChrIdx) {
			isRegionBeyond = false;
		} else {
			isRegionBeyond = (region.getStart() > read.getAlignmentEnd());
		}
		
		return isRegionBeyond;
	}
}
