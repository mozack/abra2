/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.HashMap;
import java.util.Map;

import net.sf.samtools.MyReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

public class PairValidatingRealignmentWriter implements RealignmentWriter {

	private SAMFileWriter writer;
	private ReAligner realigner;
	
	// key = read name, value = Reads
	private Map<String, Reads> firstInPair = new HashMap<String, Reads>();
	private Map<String, Reads> secondInPair = new HashMap<String, Reads>();
	
	private int realignCount = 0;
	
	private static final int INSERT_THRESHOLD = 5000;
	
	public PairValidatingRealignmentWriter(ReAligner realigner, SAMFileWriter writer) {
		
		this.writer = writer;
		this.realigner = realigner;
		throw new UnsupportedOperationException("don't run this...");
	}
	
	static long count = 1;
	
	private void put(Map<String, Reads> map, Reads reads) {
		String key = reads.getOrigRead().getReadName();
//		reads.stringify();
		map.put(key, reads);
	}
	
	public void addAlignment(SAMRecord updatedRead, SAMRecord origRead) {
		
//		if (updatedRead != null) {
//			updatedRead.clearAttributes();
//		}
//		
//		origRead.clearAttributes();
		
		if (updatedRead == null) {
			output(new Reads(updatedRead, origRead));
		} else if ((Math.abs(origRead.getInferredInsertSize()) > INSERT_THRESHOLD) || (origRead.getInferredInsertSize() == 0)) {
			output(new Reads(updatedRead, origRead));
		} else if (origRead.getFirstOfPairFlag()) {
			Reads first = new Reads(updatedRead, origRead);
			Reads second = secondInPair.get(origRead.getReadName());
			if (second != null) {
				outputPair(first, second);
				secondInPair.remove(origRead.getReadName());
			} else {
				put(firstInPair, first);
				//firstInPair.put(origRead.getReadName(), first);
			}
		} else if (origRead.getSecondOfPairFlag()) {
			Reads first = firstInPair.get(origRead.getReadName());
			Reads second = new Reads(updatedRead, origRead);
			if (first != null) {
				outputPair(first, second);
				firstInPair.remove(origRead.getReadName());
			} else {
				put(secondInPair, second);
				//secondInPair.put(origRead.getReadName(), second);
			}
		} else {
			// Unpaired read.  Just output it.
			output(new Reads(updatedRead, origRead));
		}
		
		if ((count++ % 100000) == 0) {
			System.out.println("firstInPair size: " + firstInPair.size());
			System.out.println("secondInPair size: " + secondInPair.size());
		}
		/*
		
		// If insert size exceed the threshold, just output what we have.
		// Use the updated read if available.  if not, just output the original read
		if ((Math.abs(origRead.getInferredInsertSize()) > INSERT_THRESHOLD) || (origRead.getInferredInsertSize() == 0)) {
			if (updatedRead != null) {
				writer.addAlignment(updatedRead);
			} else {
				realigner.adjustForStrand(contigAlignedRead.getReadNegativeStrandFlag(), origRead);
				writer.addAlignment(origRead);
			}
		} else {
			if (origRead.getFirstOfPairFlag()) {
				Reads first = new Reads(contigAlignedRead, updatedRead, origRead);
				Reads second = secondInPair.get(origRead.getReadName());
				if (second != null) {
					outputPair(first, second);
					secondInPair.remove(origRead.getReadName());
				} else {
					firstInPair.put(origRead.getReadName(), first);
				}
			} else if (origRead.getSecondOfPairFlag()) {
				Reads first = firstInPair.get(origRead.getReadName());
				Reads second = new Reads(contigAlignedRead, updatedRead, origRead);
				if (first != null) {
					outputPair(first, second);
					firstInPair.remove(origRead.getReadName());
				} else {
					secondInPair.put(origRead.getReadName(), second);
				}
			} else {
				// Unpaired read.  Just output it.
				output(new Reads(contigAlignedRead, updatedRead, origRead));
			}
		}
		*/
	}
	
	private void outputPair(Reads first, Reads second) {
		checkInsertLength(first, second);
		output(first);
		output(second);
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
	
	private void checkInsertLength(Reads first, Reads second) {
		
		boolean isDone = false;
		
		if ((first.getUpdatedRead() != null) && (second.getUpdatedRead() != null)) {
			if (isSameChromosome(first.getUpdatedRead(), second.getUpdatedRead())) {
				int len = getInsertLength(first.getUpdatedRead(), second.getUpdatedRead());
				
				if (len < 2 * INSERT_THRESHOLD) {
					isDone = true;
				}
			}
		}
		
		if (!isDone) {
			if (first.getUpdatedRead() != null) {
				// Compare first updated read to second orig read
				int len = getInsertLength(first.getUpdatedRead(), second.getOrigRead());
				if (len < 2 * INSERT_THRESHOLD) {
					second.clearUpdatedRead();
					isDone = true;
				}
			}
		}
		
		if (!isDone) {
			if (second.getUpdatedRead() != null) {
				// Compare second updated read to first orig read
				int len = getInsertLength(first.getOrigRead(), second.getUpdatedRead());
				if (len < 2 * INSERT_THRESHOLD) {
					first.clearUpdatedRead();
					isDone = true;
				}
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
	
	private void output(Reads reads) {
		if (reads.getUpdatedRead() != null) {
			if (reads.getUpdatedRead().getAttribute("YO") != null) {
				realignCount += 1;
			}
			writer.addAlignment(reads.getUpdatedRead());
		} else {
			SAMRecord orig = reads.getOrigRead();
//			realigner.adjustForStrand(reads.getContigAlignedRead().getReadNegativeStrandFlag(), orig);
			writer.addAlignment(orig);
		}
	}
	
	public int flush() {
		// Mate did not realign. Use original alignment info to check insert length
		for (Reads reads : firstInPair.values()) {
			checkOrigAndOutput(reads);
		}
		
		for (Reads reads : secondInPair.values()) {
			checkOrigAndOutput(reads);
		}
		
		return realignCount;
	}
	
	private void checkOrigAndOutput(Reads reads) {
		if (!reads.getUpdatedRead().getReferenceName().equals(reads.getOrigRead().getReferenceName())) {
			reads.clearUpdatedRead();
		} else if (Math.abs(reads.getUpdatedRead().getAlignmentStart() - reads.getOrigRead().getAlignmentStart()) > 2 * INSERT_THRESHOLD) {
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
				updatedRead = MyReader.getRead(updatedReadStr, realigner.getHeader());
			}
			
			return updatedRead;
		}

		public SAMRecord getOrigRead() {
			if (origRead == null) {
				origRead = MyReader.getRead(origReadStr, realigner.getHeader());
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
}
