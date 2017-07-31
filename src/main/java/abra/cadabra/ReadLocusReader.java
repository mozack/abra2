package abra.cadabra;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import abra.Feature;
import abra.SAMRecordUtils;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

public class ReadLocusReader implements Iterable<ReadsAtLocus> {

	private SamReader samReader;
	private Feature region;
	
	public ReadLocusReader(String samFile) {
		this(samFile, null);
	}
	
	public ReadLocusReader(String samFile, Feature region) {
        samReader = SAMRecordUtils.getSamReader(samFile);
        this.region = region;
	}
	
	@Override
	public Iterator<ReadsAtLocus> iterator() {
		return new ReadLocusIterator(samReader, region);
	}
	
	public SAMFileHeader getSamHeader() {
		return samReader.getFileHeader(); 
	}

	private static class ReadLocusIterator implements Iterator<ReadsAtLocus> {
		
		private Iterator<SAMRecord> samIter;
		private String currentChr = "";
		private int currentPos = -1;
		private List<SAMRecord> readCache = new ArrayList<SAMRecord>();
		private ReadsAtLocus nextCache;
		
		public ReadLocusIterator(SamReader samReader, Feature region) {
	        	  
			if (region != null) {
				samIter = new ForwardShiftInsertIterator(samReader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd()));
			} else {
				samIter = new ForwardShiftInsertIterator(samReader.iterator());
			}
		}
		
		@Override
		public boolean hasNext() {
			nextCache = next();
			return nextCache != null;
		}

		@Override
		public ReadsAtLocus next() {
			
			if (nextCache != null) {
				// Return the cached value
				ReadsAtLocus ret = nextCache;
				nextCache = null;
				return ret;
			}
			
			List<SAMRecord> reads = new ArrayList<SAMRecord>();
			
			loadReadsIntoCache();
			
			boolean isLocusAdvanced = getCachedReadsAtCurrentLocus(reads);
			
			if (isLocusAdvanced) {
				getCachedReadsAtCurrentLocus(reads);
			}
			ReadsAtLocus readsAtLocus = null;
			if (!reads.isEmpty()) {
				readsAtLocus = new ReadsAtLocus(currentChr, currentPos, reads);
			}
			
			currentPos += 1;
									
			return readsAtLocus;
		}
		
		private void loadReadsIntoCache() {
			boolean shouldReadFromFile = false;
			
			if (readCache.isEmpty()) {
				shouldReadFromFile = true;
			}
			else {
				SAMRecord last = readCache.get(readCache.size()-1);
				if (getAlignmentStart(last) <= currentPos && last.getReferenceName().equals(currentChr)) {
					shouldReadFromFile = true;
				}
			}
			
			while (shouldReadFromFile && samIter.hasNext()) {
				SAMRecord read = samIter.next();

				// Skip over unmapped reads
				if (!read.getReadUnmappedFlag()) {
				
					if (readCache.isEmpty() && !read.getReferenceName().equals(currentChr)) {
						currentChr = read.getReferenceName();
						currentPos = getAlignmentStart(read);
					}
					
					readCache.add(read);
					
					if (getAlignmentStart(read) > currentPos || !read.getReferenceName().equals(currentChr)) {
						shouldReadFromFile = false;
					}
				}
			}			
		}
		
		private int getAlignmentStart(SAMRecord read) {
			int start = read.getAlignmentStart();
			
			if (read.getCigarLength() > 0 && read.getCigar().getCigarElement(0).getOperator() == CigarOperator.I) {
				start = start - 1;
			}
			
			return start;
		}
		
		// Returns true if current position is advanced to new locus
		private boolean getCachedReadsAtCurrentLocus(List<SAMRecord> reads) {
			
			reads.clear();
			
			Iterator<SAMRecord> cacheIter = readCache.iterator();
			
			String nextChr = null;
			int nextPos = -1;
			
			while (cacheIter.hasNext()) {
				SAMRecord read = cacheIter.next();
				
				if (read.getAlignmentEnd() < currentPos && read.getReferenceName().equals(currentChr)) {
					// We've gone past the end of this read, so remove from cache.
					cacheIter.remove();
				} else if (getAlignmentStart(read) <= currentPos && read.getAlignmentEnd() >= currentPos) {
					// This read spans the current locus of interest.
					reads.add(read);
				} else {
					// This read is beyond the current locus.
					if (nextChr == null) {
						nextChr = read.getReferenceName();
						nextPos = getAlignmentStart(read);
					}
				}
			}
			
			if (reads.isEmpty() && nextChr != null) {
				currentChr = nextChr;
				currentPos = nextPos;
				
				return true;
			}
			
			return false;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
	
	public static void main(String[] args) {
		String file = "/home/lmose/dev/abra/dream/small.sort.bam";
		
		ReadLocusReader r = new ReadLocusReader(file, null);
		
		for (ReadsAtLocus readsAtLocus : r) {
			System.out.println(readsAtLocus);
		}
	}
}
