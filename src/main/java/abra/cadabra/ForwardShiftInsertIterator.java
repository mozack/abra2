package abra.cadabra;

import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ForwardShiftInsertIterator implements Iterator<SAMRecord> {
	
	private Iterator<SAMRecord> iter;
	
	// Cached reads sorted by adjusted alignment start position.
	private NavigableSet<InsertShiftSAMRecord> cache = new TreeSet<InsertShiftSAMRecord>();

	public ForwardShiftInsertIterator(Iterator<SAMRecord> iter) {
		this.iter = iter;
	}
	
	@Override
	public boolean hasNext() {
		return !cache.isEmpty() || iter.hasNext();
	}

	@Override
	public SAMRecord next() {

		SAMRecord read = null;
		boolean isCacheUpToDate = false;
		if (!cache.isEmpty()) {
			InsertShiftSAMRecord first = cache.first();
			InsertShiftSAMRecord last = cache.last();
			
			// Don't seek too far ahead
			if (last.read.getAlignmentStart() > first.read.getAlignmentStart()+2 || !last.read.getReferenceName().equals(first.read.getReferenceName())) {
				isCacheUpToDate = true;
			}
			
			read = first.read;
		} else {
			read = iter.next();
			cache.add(new InsertShiftSAMRecord(read));
		}
		
		int cacheStart = read.getAlignmentStart() + 1;
		String cacheChromosome = read.getReferenceName();
		
		while (!isCacheUpToDate && iter.hasNext() && read.getAlignmentStart() <= cacheStart+1 && read.getReferenceName().equals(cacheChromosome)) {
			read = iter.next();
			cache.add(new InsertShiftSAMRecord(read));
		}
		
		return cache.pollFirst().read;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	
	static class InsertShiftSAMRecord implements Comparable<InsertShiftSAMRecord> {
		private SAMRecord read;
		
		InsertShiftSAMRecord(SAMRecord read) {
			this.read = read;
		}
		
		public int getAlignmentStart() {
			int start = read.getAlignmentStart();
			
			if (read.getCigarLength() > 0 && read.getCigar().getCigarElement(0).getOperator() == CigarOperator.I) {
				start = start - 1;
			}
			
			return start;
		}

		@Override
		public int compareTo(InsertShiftSAMRecord that) {
			int compare = this.read.getReferenceIndex() - that.read.getReferenceIndex();
			if (compare == 0) {
				compare = this.getAlignmentStart() - that.getAlignmentStart();
			}
			if (compare == 0) {
				compare = this.read.getReadName().compareTo(that.read.getReadName());
			}
			if (compare == 0) {
				compare = this.read.getFlags() - that.read.getFlags();
			}
			return compare;
		}
	}
}
