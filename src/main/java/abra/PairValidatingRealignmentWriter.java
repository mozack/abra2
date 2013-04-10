package abra;

import java.util.Map;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

public class PairValidatingRealignmentWriter implements RealignmentWriter {

	private SAMFileWriter writer;
	private ReAligner realigner;
	
	// key = read name, value = Reads
	private Map<String, Reads> firstInPair;
	private Map<String, Reads> secondInPair;
	
	private int realignCount = 0;
	
	private static final int INSERT_THRESHOLD = 5000;
	
	public PairValidatingRealignmentWriter(ReAligner realigner, SAMFileWriter writer) {
		this.writer = writer;
		this.realigner = realigner;
	}
	
	public void addAlignment(SAMRecord contigAlignedRead, SAMRecord updatedRead, SAMRecord origRead) {
		
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
	}
	
	private void outputPair(Reads first, Reads second) {
		checkInsertLength(first, second);
		output(first);
		output(second);
	}
	
	private int getInsertLength(SAMRecord read1, SAMRecord read2) {
		int start = Math.min(read1.getAlignmentStart(), read2.getAlignmentStart());
		int end = Math.max(read1.getAlignmentEnd(), read2.getAlignmentEnd());
		int len = end - start;
		
		return len;
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
			realigner.adjustForStrand(reads.getContigAlignedRead().getReadNegativeStrandFlag(), orig);
			writer.addAlignment(orig);
		}
	}
	
	public int flush() {
		for (Reads reads : firstInPair.values()) {
			output(reads);
		}
		
		for (Reads reads : secondInPair.values()) {
			output(reads);
		}
		
		return realignCount;
	}
	
	static class Reads {
		private SAMRecord contigAlignedRead;
		private SAMRecord updatedRead;
		private SAMRecord origRead;
		
		public Reads(SAMRecord contigAlignedRead, SAMRecord updatedRead,
				SAMRecord origRead) {
			this.contigAlignedRead = contigAlignedRead;
			this.updatedRead = updatedRead;
			this.origRead = origRead;
		}

		public SAMRecord getContigAlignedRead() {
			return contigAlignedRead;
		}

		public SAMRecord getUpdatedRead() {
			return updatedRead;
		}

		public SAMRecord getOrigRead() {
			return origRead;
		}
		
		public void clearUpdatedRead() {
			this.updatedRead = null;
		}
	}
}
