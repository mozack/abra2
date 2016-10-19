package abra;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class AltContigGenerator {
	
	// TODO: Multiple arbitrary cutoffs here, need to optimize
	// Phred 13 is ~5% error rate
	private static int SOFT_CLIPPING_MIN_BASE_QUAL = 13;
	
	private boolean hasHighQualitySoftClipping(SAMRecord read, int start, int length) {
		
		int numHighQualBases = 0;
		int requiredHighQualBases = (int) (.8 * length);
		
		for (int bq = start; bq < start+length; bq++) {
			if (read.getBaseQualities()[bq] >= SOFT_CLIPPING_MIN_BASE_QUAL) {
				numHighQualBases += 1;
				
				if (numHighQualBases >= requiredHighQualBases) {
					return true;
				}
			}
		}
		
		return false;
	}
	
	// Return true if read contains a soft clip element >= 5 bases long with 80% or more of bases exceeding min base qual
	private boolean hasHighQualitySoftClipping(SAMRecord read) {
		
		boolean hasHighQualitySoftClipping = false;
		
		if (read.getCigarLength() > 1) {
			// Check first cigar element
			CigarElement elem = read.getCigar().getCigarElement(0);
			if (elem.getOperator() == CigarOperator.S) {
				int elemStart = 0;
				hasHighQualitySoftClipping = hasHighQualitySoftClipping(read, elemStart, elem.getLength());
			}
			
			// Check last Cigar element
			if (!hasHighQualitySoftClipping) {
				elem = read.getCigar().getCigarElement(read.getCigarLength()-1);
				if (elem.getOperator() == CigarOperator.S) {
					int elemStart = read.getReadLength() - elem.getLength();
					hasHighQualitySoftClipping = hasHighQualitySoftClipping(read, elemStart, elem.getLength());
				}
			}
		}
		
		return hasHighQualitySoftClipping;
	}

	public Collection<String> getAltContigs(List<List<SAMRecordWrapper>> readsList, CompareToReference2 c2r, int readLength) {
		
		Set<String> contigs = new HashSet<String>();
		
		HashSet<Indel> indels = new HashSet<Indel>();
		
		for (List<SAMRecordWrapper> reads : readsList) {
			for (SAMRecordWrapper readWrapper : reads) {
				SAMRecord read = readWrapper.getSamRecord();
				
				// TODO: Any other filters here.  Higher MQ threshold?
				if (read.getMappingQuality() > 0) {
					
					// For now only use indels bracketed by 2 M elements
					// TODO: Handle clipping / complex indels
					List<CigarElement> elems = read.getCigar().getCigarElements();
					if (elems.size() == 3 && 
						elems.get(0).getOperator() == CigarOperator.M && 
						elems.get(2).getOperator() == CigarOperator.M &&
						(elems.get(1).getOperator() == CigarOperator.D || elems.get(1).getOperator() == CigarOperator.I)) {
						
						String insertBases = null;
						char type = '0';
						if (elems.get(1).getOperator() == CigarOperator.D) {
							type = 'D';
						} else if (elems.get(1).getOperator() == CigarOperator.I) {
							type = 'I';
	//						System.out.println("read: " + read.getSAMString() + ", elem0: " + elems.get(0).getLength() + ", elem1: " + elems.get(1).getLength());
							int start = elems.get(0).getLength();
							int stop =  start + elems.get(1).getLength();
							insertBases = read.getReadString().substring(start, stop);
						}
						
						Indel indel = new Indel(type, read.getReferenceName(), read.getAlignmentStart() + elems.get(0).getLength(), elems.get(1).getLength(), insertBases);
						indels.add(indel);
					}
					
					// Add high quality soft clipped reads
					if (hasHighQualitySoftClipping(readWrapper.getSamRecord())) {
						System.err.println("SOFT_SEQ: " + readWrapper.getSamRecord());
						contigs.add(read.getReadString());
					}
				}
			}
		}
		
		// Current set of contigs is from soft clipping.  These can grow to be too big.  Downsample if necessary.
		// TODO: Is this really the right thing to do?
		int maxSoftClippedContigs = 100;
		if (contigs.size() > maxSoftClippedContigs) {
			List<String> contigList = new ArrayList<String>(contigs);
			Collections.shuffle(contigList);
			contigs = new HashSet<String>(contigList.subList(0, maxSoftClippedContigs));
		}
		
		for (Indel indel : indels) {
			if (indel.type == 'D') {
				// Pull in read length sequence from both sides of deletion.
				int leftStart = indel.pos - readLength;
				int rightStart = indel.pos + indel.length;
				String leftSeq = c2r.getSequence(indel.chr, leftStart, readLength);
				String rightSeq = c2r.getSequence(indel.chr, rightStart, readLength);
				String seq = leftSeq + rightSeq;
				System.err.println("INDEL_SEQ: chr: " + indel.chr + "  pos: " + indel.pos + " type: " + indel.type + " length: " + indel.length + " bases: " + indel.insert + " seq: [" + seq + "]");
				
				contigs.add(seq);
			} else if (indel.type == 'I') {
				// Pull in read length sequence on both sides of insertion
				int leftStart = indel.pos - readLength;
				int rightStart = indel.pos;
				String leftSeq = c2r.getSequence(indel.chr, leftStart, readLength);
				String rightSeq = c2r.getSequence(indel.chr, rightStart, readLength);
				String seq = leftSeq + indel.insert + rightSeq;
				System.err.println("INDEL_SEQ: chr: " + indel.chr + "  pos: " + indel.pos + " type: " + indel.type + " length: " + indel.length + " bases: " + indel.insert + " seq: [" + seq + "]");
				
				contigs.add(seq);
			}
		}
		
		return contigs;
	}
	
	static class Indel {
		char type;
		String chr;
		int pos;
		int length;
		String insert;
		
		Indel(char type, String chr, int pos, int length, String insert) {
			this.type = type;
			this.chr = chr;
			this.pos = pos;
			this.length = length;
			this.insert = insert;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((chr == null) ? 0 : chr.hashCode());
			result = prime * result
					+ ((insert == null) ? 0 : insert.hashCode());
			result = prime * result + length;
			result = prime * result + pos;
			result = prime * result + type;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Indel other = (Indel) obj;
			if (chr == null) {
				if (other.chr != null)
					return false;
			} else if (!chr.equals(other.chr))
				return false;
			if (insert == null) {
				if (other.insert != null)
					return false;
			} else if (!insert.equals(other.insert))
				return false;
			if (length != other.length)
				return false;
			if (pos != other.pos)
				return false;
			if (type != other.type)
				return false;
			return true;
		}
	}
}
