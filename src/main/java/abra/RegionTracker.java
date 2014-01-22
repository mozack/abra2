package abra;

import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * Keeps track of current region and provides
 * utility methods on reads related to that region.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class RegionTracker {
	
	private Feature currentRegion;
	private Iterator<Feature> regionIter;
	private SAMFileHeader header;

	public RegionTracker(List<Feature> regions, SAMFileHeader header) {
		regionIter = regions.iterator();
        if (regionIter.hasNext()) {
        	currentRegion = regionIter.next();
        }
		this.header = header;
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
