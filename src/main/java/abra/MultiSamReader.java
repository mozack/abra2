package abra;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

public class MultiSamReader implements Iterable<SAMRecordWrapper> {
	
	//TODO: Upgrade to newer implementation
	private SamReader[] readers;
	private Iterator<SAMRecord>[] iterators;
	private SAMRecordWrapper[] nextRecord;
	private int minMapqForAssembly;
	private boolean isPairedEnd;
	
	// Iterator used by clients
	private Iterator<SAMRecordWrapper> clientIterator;
	
	public MultiSamReader(String[] inputBams, int minMapqForAssembly, boolean isPairedEnd, Feature region) {
		
		//TODO: Assert all SAM Headers have same sequence dict
		readers = new SamReader[inputBams.length];
		nextRecord = new SAMRecordWrapper[inputBams.length];
		iterators = new Iterator[inputBams.length];
		this.minMapqForAssembly = minMapqForAssembly;
		this.isPairedEnd = isPairedEnd;
		
		int idx = 0;
		for (String bamFileName : inputBams) {
			SamReader reader = SAMRecordUtils.getSamReader(bamFileName);
			
			readers[idx] = reader;
			
			// TODO: Pad by region size?
			iterators[idx] = readers[idx].queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
			
			// cache next record
			cacheNextRecord(idx);
			
			idx += 1;
		}
		
		clientIterator = new MultiSamReaderIterator(this);
	}
	
	private void cacheNextRecord(int sampleIdx) {
		nextRecord[sampleIdx] = getNext(sampleIdx);
	}
	
	public SAMFileHeader getSAMFileHeader() {
		return readers[0].getFileHeader();
	}
	
	public void close() throws IOException {
		for (SamReader reader : readers) {
			reader.close();
		}
	}
	
	private boolean isFiltered(SAMRecord read) {
		return SAMRecordUtils.isFiltered(isPairedEnd, read);
	}
	
	private SAMRecordWrapper getNext(int idx) {
		SAMRecordWrapper record = null;
		if (iterators[idx].hasNext()) {
			SAMRecord read = iterators[idx].next();
			// If no genomic location is assigned, we've reached the unmapped read pairs.  Do not continue...
			// TODO: Need to include these in final bam files
			if (read.getReferenceIndex() >= 0) {
				record = new SAMRecordWrapper(read, isFiltered(read), shouldAssemble(read), idx);
			}
		}
		
		return record;
	}
	
	private boolean shouldAssemble(SAMRecord read) {
		return (
			(!read.getDuplicateReadFlag()) && 
			(!read.getReadFailsVendorQualityCheckFlag()) &&
			(read.getMappingQuality() >= this.minMapqForAssembly || read.getReadUnmappedFlag()) &&
			SAMRecordUtils.isPrimary(read));  // Was previously an id check, so supplemental / secondary alignments could be included
	}
	
	@Override
	public Iterator<SAMRecordWrapper> iterator() {
		return clientIterator;
	}
	
	static class MultiSamReaderIterator implements Iterator<SAMRecordWrapper> {
		
		private MultiSamReader multiSamReader;
		
		MultiSamReaderIterator(MultiSamReader multiSamReader) {
			this.multiSamReader = multiSamReader;
		}

		@Override
		public boolean hasNext() {
			// Return true if any sample has another read
			for (SAMRecordWrapper record : multiSamReader.nextRecord) {
				if (record != null) {
					return true;
				}
			}
			
			return false;
		}

		@Override
		public SAMRecordWrapper next() {
			// Return the first read across samples by genomic coordinate
			SAMRecordWrapper nextRecord = null;
			int bestChr = Integer.MAX_VALUE;
			int bestPos = Integer.MAX_VALUE;
			int bestSampleIdx = -1;
			
			for (int i=0; i<multiSamReader.nextRecord.length; i++) {
				
				SAMRecordWrapper record = multiSamReader.nextRecord[i];
				
				if (record != null) {
					if ( (record.getSamRecord().getReferenceIndex() < bestChr) || 
						 (record.getSamRecord().getReferenceIndex() == bestChr && record.getSamRecord().getAlignmentStart() < bestPos)) {
						nextRecord = record;
						bestChr = record.getSamRecord().getReferenceIndex();
						bestPos = record.getSamRecord().getAlignmentStart();
						bestSampleIdx = i;
					}
				}
			}

			// Replace current read in cache
			multiSamReader.cacheNextRecord(bestSampleIdx);
			
			return nextRecord;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
}
