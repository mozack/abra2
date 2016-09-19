package abra;

import java.io.File;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

public class MultiSamReader implements Iterable<SAMRecordWrapper> {
	
	//TODO: Upgrade to newer implementation
	private SAMFileReader[] readers;
	private Iterator<SAMRecord>[] iterators;
	private SAMRecordWrapper[] nextRecord;
	private ReAligner realigner;
	
	// Iterator used by clients
	private Iterator<SAMRecordWrapper> clientIterator;

	public MultiSamReader(String[] inputBams, ReAligner realigner) {
		
		//TODO: Assert all SAM Headers have same sequence dict
		readers = new SAMFileReader[inputBams.length];
		nextRecord = new SAMRecordWrapper[inputBams.length];
		iterators = new Iterator[inputBams.length];
		this.realigner = realigner;
		
		int idx = 0;
		for (String bamFileName : inputBams) {
			SAMFileReader reader = new SAMFileReader(new File(bamFileName));
			reader.setValidationStringency(ValidationStringency.SILENT);
			
			readers[idx] = reader;
			iterators[idx] = readers[idx].iterator();
			
			// cache next record
			nextRecord[idx] = getNext(idx);
			
			idx += 1;
		}
		
		clientIterator = new MultiSamReaderIterator(this);
	}
	
	public SAMFileHeader getSAMFileHeader() {
		return readers[0].getFileHeader();
	}
	
	public void close() {
		for (SAMFileReader reader : readers) {
			reader.close();
		}
	}
	
	private SAMRecordWrapper getNext(int idx) {
		SAMRecordWrapper record = null;
		if (iterators[idx].hasNext()) {
			SAMRecord read = iterators[idx].next();
			// If no genomic location is assigned, we've reached the unmapped read pairs.  Do not continue...
			if (read.getReferenceIndex() >= 0) {
				record = new SAMRecordWrapper(read, realigner.isFiltered(read), shouldAssemble(read), idx);
			}
		}
		
		return record;
	}
	
	private boolean shouldAssemble(SAMRecord read) {
		return ((!read.getDuplicateReadFlag()) && 
			(!read.getReadFailsVendorQualityCheckFlag()) &&
			(read.getMappingQuality() >= realigner.getMinMappingQuality() || read.getReadUnmappedFlag()) &&
			(!read.getNotPrimaryAlignmentFlag()));  // Was previously an id check, so supplemental / secondary alignments could be included
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
			
			for (SAMRecordWrapper record : multiSamReader.nextRecord) {
				if ( (record.getSamRecord().getReferenceIndex() < bestChr) || 
					 (record.getSamRecord().getReferenceIndex() == bestChr && record.getSamRecord().getAlignmentStart() == bestPos)) {
					nextRecord = record;
					bestChr = record.getSamRecord().getReferenceIndex();
					bestPos = record.getSamRecord().getAlignmentStart();
				}
			}
			
			return nextRecord;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
}
