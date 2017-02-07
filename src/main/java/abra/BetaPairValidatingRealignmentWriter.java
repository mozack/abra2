/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMLineParser;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;

/**
 * Manages writing paired reads to final output.  Reads that were previously in
 * a "proper pair" will not be modified to become an "improper pair".
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class BetaPairValidatingRealignmentWriter implements RealignmentWriter {

	private SAMFileWriter writer;
//	private ReAligner realigner;
	private IndelShifter indelShifter = new IndelShifter();
	
	private int realignCount = 0;
	
	private int maxInsertLength;
	private int minInsertLength;
	
	private String candidatesSam;
	private SAMFileWriter candidatesSamWriter;
	private SAMLineParser parser;
	
	private CompareToReference2 c2r;
	private SAMFileHeader header;
	
	public BetaPairValidatingRealignmentWriter(CompareToReference2 c2r,
			SAMFileWriter writer, String tempDir, int minInsertLen, int maxInsertLen) {
		this.writer = writer;
		this.c2r = c2r;
		
		header = writer.getFileHeader().clone();
		header.setSortOrder(SortOrder.queryname);
		
		parser = new SAMLineParser(new DefaultSAMRecordFactory(),
                ValidationStringency.SILENT, header,
                null, null);
		
		candidatesSam = tempDir + "/candidates.bam";
		
		candidatesSamWriter = new SAMFileWriterFactory().makeBAMWriter(
				header, false, new File(candidatesSam), 1);
		
		this.minInsertLength = minInsertLen;
		this.maxInsertLength = maxInsertLen;
	}
	
	long count = 1;
	int numCandidates = 0;

	private boolean isValidInsertLength(int insertLen) {
		return Math.abs(insertLen) >= minInsertLength && insertLen <= maxInsertLength;
	}
	
	private boolean isValidOrientation(SAMRecord read1, int read2Start, boolean isRead2OnReverseStrand) {
		boolean isFirstReadOnReverseStrand;
		boolean isSecondReadOnReverseStrand;
		
		if (read1.getAlignmentStart() < read2Start) {
			isFirstReadOnReverseStrand = read1.getReadNegativeStrandFlag();
			isSecondReadOnReverseStrand = isRead2OnReverseStrand;
		} else {
			isFirstReadOnReverseStrand = isRead2OnReverseStrand;
			isSecondReadOnReverseStrand = read1.getReadNegativeStrandFlag();
		}
				
		return !isFirstReadOnReverseStrand && isSecondReadOnReverseStrand;
	}
	
	private boolean isValidOrientation(SAMRecord read1, SAMRecord read2) {
		return isValidOrientation(read1, read2.getAlignmentStart(), read2.getReadNegativeStrandFlag());
	}
	
	public void addAlignment(SAMRecord updatedRead, SAMRecord origRead) {
		
		
		if (updatedRead == null) {
			// Just output the original read
			output(new Reads(updatedRead, origRead));
		} else if (updatedRead.getAttribute("YO") == null) {
			// Updated read has not moved, just output it
			output(new Reads(updatedRead, origRead));
		} else if ((!origRead.getReadPairedFlag()) || (!origRead.getProperPairFlag())) {
			// Original read not part of "proper pair", output updated read
			output(new Reads(updatedRead, origRead));
		} else {
			// Output candidate to temp bam for comparison with mate
			writeToTempFile(candidatesSamWriter, updatedRead, origRead);
			numCandidates += 1;
		}
		
		if ((count++ % 100000) == 0) {
			System.err.println("Num candidates: " + numCandidates);
		}
	}
	
	private void writeToTempFile(SAMFileWriter writer, SAMRecord updatedRead, SAMRecord origRead) {
		updatedRead.setAttribute("YG", origRead.getSAMString());
		writer.addAlignment(updatedRead);
	}
	
	private void outputPair(Reads first, Reads second) {
		checkPairValidity(first, second);
		
		if ((first.getUpdatedRead() != null) && (second.getUpdatedRead() != null)) {
			
			int insertLen = getInsertLength(first.getUpdatedRead(), second.getUpdatedRead());
			
			// Both reads are realigned.  insert length and read orientation is proper.
			first.getUpdatedRead().setProperPairFlag(true);
			first.getUpdatedRead().setMateUnmappedFlag(false);
			first.getUpdatedRead().setMateAlignmentStart(second.getUpdatedRead().getAlignmentStart());
			first.getUpdatedRead().setMateReferenceName(second.getUpdatedRead().getReferenceName());
			first.getUpdatedRead().setMateNegativeStrandFlag(second.getUpdatedRead().getReadNegativeStrandFlag());
			
			second.getUpdatedRead().setProperPairFlag(true);
			second.getUpdatedRead().setMateUnmappedFlag(false);
			second.getUpdatedRead().setMateAlignmentStart(first.getUpdatedRead().getAlignmentStart());
			second.getUpdatedRead().setMateReferenceName(first.getUpdatedRead().getReferenceName());
			second.getUpdatedRead().setMateNegativeStrandFlag(first.getUpdatedRead().getReadNegativeStrandFlag());
			
			if (first.getUpdatedRead().getAlignmentStart() < second.getUpdatedRead().getAlignmentStart()) {
				first.getUpdatedRead().setInferredInsertSize(insertLen);
				second.getUpdatedRead().setInferredInsertSize(-insertLen);
			} else {
				first.getUpdatedRead().setInferredInsertSize(-insertLen);
				second.getUpdatedRead().setInferredInsertSize(insertLen);				
			}
		}
		
		output(first);
		output(second);
	}
	
	private int getInsertGap(int read1Start, int read1End, int read2Start, int read2End) {
		int gap = 0;
		
		if (read1Start < read2Start) {
			gap = read2Start - read1End;
		} else {
			gap = read1Start - read2End;
		}
		
		return gap;
	}
	
	private int getInsertLength(int read1Start, int read1End, int read2Start, int read2End) {
		int start = Math.min(read1Start, read2Start);
		int end = Math.max(read1End, read2End);
		int len = end - start;
		
		return len;
	}
	
	private int getInsertLength(SAMRecord read1, SAMRecord read2) {
		return getInsertLength(read1.getAlignmentStart(), read1.getAlignmentEnd(),
				read2.getAlignmentStart(), read2.getAlignmentEnd());
	}
	
	private boolean isPairValid(SAMRecord read1, SAMRecord read2) {
		boolean isValid = false;
		
		if ((read1 != null) && (read2 != null)) {
			if (isSameChromosome(read1, read2)) {
				//int len = getInsertLength(read1, read2);
				int len = getInsertGap(read1.getAlignmentStart(), read1.getAlignmentEnd(), read2.getAlignmentStart(),
						read2.getAlignmentEnd()) +
						read1.getReadLength() + read2.getReadLength();
								
				isValid = (isValidInsertLength(len)) && (isValidOrientation(read1, read2));
			}
		}
		
		return isValid;
	}
	
	private boolean isPairValidWithOriginalMate(SAMRecord read) {
		boolean isValid = false;
		
		String mateChr = read.getMateReferenceName();
		int mateStart = read.getMateAlignmentStart();
		int mateEnd = mateStart + read.getReadLength();
		boolean isMateOnReverseStrand = read.getMateNegativeStrandFlag();
		
		if (read.getReferenceName().equals(mateChr)) {
//			int len = getInsertLength(read.getAlignmentStart(), read.getAlignmentEnd(),
//					mateStart, mateEnd);
			
			//TODO: mateEnd may be slightly off here.
			int len = getInsertGap(read.getAlignmentStart(), read.getAlignmentEnd(), mateStart,
					mateEnd) + read.getReadLength() + read.getReadLength();
			
			isValid = isValidInsertLength(len) && isValidOrientation(read, mateStart, isMateOnReverseStrand);
		}
		
		return isValid;
	}
	
	private void checkPairValidity(Reads first, Reads second) {
		
		boolean isDone = false;
		
		if (isPairValid(first.getUpdatedRead(), second.getUpdatedRead())) {
			isDone = true;
		}
				
		if (!isDone) {
			if (isPairValid(first.getUpdatedRead(), second.getOrigRead())) {
				second.clearUpdatedRead();
				isDone = true;				
			}			
		}
		
		if (!isDone) {
			if (isPairValid(first.getOrigRead(), second.getUpdatedRead())) {
				first.clearUpdatedRead();
				isDone = true;
			}
		}
		
		if (!isDone) {
			first.clearUpdatedRead();
			second.clearUpdatedRead();
		}
	}
	
	private boolean isSameChromosome(SAMRecord read1, SAMRecord read2) {
		return read1.getReferenceName().equals(read2.getReferenceName());
	}
	
	int updatedCount = 0;
	int origCount = 0;
	
	private void output(Reads reads) {
		if (reads.getUpdatedRead() != null) {
			if (reads.getUpdatedRead().getAttribute("YO") != null) {
				realignCount += 1;
			}
			addAlignment(reads.getUpdatedRead());
			updatedCount += 1;
		} else {
			SAMRecord orig = reads.getOrigRead();
			addAlignment(orig);
			origCount += 1;
		}
	}
	
	private void addAlignment(SAMRecord read) {
//		writer.addAlignment(indelShifter.shiftIndelsLeft(read, c2r));
		writer.addAlignment(read);
	}
	
	private void processCandidates() {
		System.err.println("Processing candidates");
		SimpleSamReadPairReader reader = new SimpleSamReadPairReader(candidatesSam);
		
		for (ReadPair pair : reader) {
			SAMRecord updatedRead1 = pair.getRead1();
			SAMRecord updatedRead2 = pair.getRead2();
			SAMRecord origRead1 = getOriginalRead(updatedRead1);
			SAMRecord origRead2 = getOriginalRead(updatedRead2);
			
			if ((updatedRead1 !=  null) && (updatedRead2 != null)) {
				Reads reads1 = new Reads(updatedRead1, origRead1);
				Reads reads2 = new Reads(updatedRead2, origRead2);
				
				outputPair(reads1, reads2);
			} else if (updatedRead1 != null) {
				Reads reads1 = new Reads(updatedRead1, origRead1);
				checkOrigAndOutput(reads1);
			} else if (updatedRead2 != null) {
				Reads reads2 = new Reads(updatedRead2, origRead2);
				checkOrigAndOutput(reads2);
			}
		}
		System.err.println("Done processing candidates");
	}
	
	private SAMRecord getOriginalRead(SAMRecord read) {
		SAMRecord orig = null;
		
		if (read != null) {
			String origStr = (String) read.getAttribute("YG");
			orig = parser.parseLine(origStr);
			read.setAttribute("YG", null);
		}
		
		return orig;
	}
	
	public int flush() {
		System.err.println("Flushing");
		candidatesSamWriter.close();
		processCandidates();
		
		System.err.println("updatedCount: " + updatedCount);
		System.err.println("origCount: " + origCount);
		return realignCount;
	}
	
	private void checkOrigAndOutput(Reads reads) {
		
		if (!isPairValidWithOriginalMate(reads.getUpdatedRead())) {
			reads.clearUpdatedRead();
		}
		
		output(reads);		
	}
	
	class Reads {
		private SAMRecord updatedRead;
		private SAMRecord origRead;
		
		private String updatedReadStr;
		private String origReadStr;
		
		public Reads(SAMRecord updatedRead, SAMRecord origRead) {
			this.updatedRead = updatedRead;
			this.origRead = origRead;
		}

		public SAMRecord getUpdatedRead() {
			if ( (updatedReadStr != null) && (updatedRead == null) ) {
				updatedRead = parser.parseLine(updatedReadStr);
			}
			
			return updatedRead;
		}

		public SAMRecord getOrigRead() {
			if (origRead == null) {
				origRead = parser.parseLine(origReadStr);
			}
			
			return origRead;
		}
		
		public void clearUpdatedRead() {
			this.updatedRead = null;
			this.updatedReadStr = null;
		}
		
		public void stringify() {
			if (updatedRead != null) {
				this.updatedReadStr = updatedRead.getSAMString();
			}
			this.origReadStr = origRead.getSAMString();
			
			this.updatedRead = null;
			this.origRead = null;
		}
	}
	
	public static void main(String[] args) throws IOException {
		SamReader reader = SAMRecordUtils.getSamReader("/home/lmose/dev/abra/1076/candidates.bam");
		SAMFileHeader header = reader.getFileHeader();
		reader.close();
		header.setSortOrder(SortOrder.unsorted);
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		writerFactory.setUseAsyncIo(false);
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				header, false, new File("/home/lmose/dev/abra/1076/test.bam"));

		
		BetaPairValidatingRealignmentWriter w = new BetaPairValidatingRealignmentWriter(null,
				writer, "foofi", 0, 200000);
		
		w.candidatesSam = "/home/lmose/dev/abra/1076/candidates.bam";
		w.processCandidates();
		writer.close();
	}

}
