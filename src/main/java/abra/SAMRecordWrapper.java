package abra;

import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SAMRecordWrapper {

	private SAMRecord samRecord;
	private boolean shouldAssemble;
	private boolean shouldFilter;
	private int sampleIdx;
	private boolean isUnalignedRc = false;
	private String mergedSeq = null;
	private String mergedQual = null;
	private int adjustedAlignmentStart = -1;
	private int adjustedAlignmentEnd = -1;
	
	private int bqSum = -1;
	
	public SAMRecordWrapper(SAMRecord record, boolean shouldFilter, boolean shouldAssemble, int sampleIdx) {
		this.samRecord = record;
		this.shouldFilter = true;
		this.shouldAssemble = shouldAssemble;
		this.sampleIdx = sampleIdx;
	}

	public SAMRecord getSamRecord() {
		return samRecord;
	}

	public boolean shouldAssemble() {
		return shouldAssemble;
	}

	public boolean shouldFilter() {
		return shouldFilter;
	}

	public int getSampleIdx() {
		return sampleIdx;
	}
	
	public boolean isUnalignedRc() {
		return isUnalignedRc;
	}

	public void setUnalignedRc(boolean isUnalignedRc) {
		this.isUnalignedRc = isUnalignedRc;
	}

	public void setShouldAssemble(boolean shouldAssemble) {
		this.shouldAssemble = false;
	}
	public int getAdjustedAlignmentStart() {
		
		int start = 0;
		
		if (adjustedAlignmentStart > -1) {
			start = adjustedAlignmentStart;
		} else {
		
			start = samRecord.getAlignmentStart();
			
			if (samRecord.getCigar().numCigarElements() > 0) {
				CigarElement elem = samRecord.getCigar().getCigarElement(0);
				if (elem.getOperator() == CigarOperator.S) {
					start -= elem.getLength();
					if (start < 1) {
						start = 1;
					}
				}
			}
		}
		
		return start;
	}
	
	public int getAdjustedAlignmentEnd() {
		int end = -1;
		
		if (adjustedAlignmentEnd > -1) {
			end = adjustedAlignmentEnd;
		} else {
		
			if (samRecord.getReadUnmappedFlag()) {
				end = samRecord.getAlignmentStart() + samRecord.getReadLength();
			} else {
				// Use standard alignment end and pad for soft clipping if necessary
				end = samRecord.getAlignmentEnd();
				
				if (samRecord.getCigar().numCigarElements() > 0) {
					CigarElement elem = samRecord.getCigar().getCigarElement(samRecord.getCigar().numCigarElements()-1);
					if (elem.getOperator() == CigarOperator.S) {
						end += elem.getLength();
					}
				}
			}
		}

		return end;
	}
	
	public int getReadLength() {
		int length = this.samRecord.getReadLength();
		if (hasMergedSeq()) {
			length = mergedSeq.length();
		}
		
		return length;
	}
	
	public String getMergedSeq() {
		return mergedSeq;
	}
	
	public String getMergedQual() {
		return mergedQual;
	}
	
	public String getSeq() {
		String seq;
		if (mergedSeq != null) {
			seq = mergedSeq;
		} else {
			seq = samRecord.getReadString();
		}
		
		return seq;
	}
	
	public String getQual() {
		String qual;
		if (mergedQual != null) {
			qual = mergedQual;
		} else {
			qual = samRecord.getBaseQualityString();
		}
		
		return qual;
	}
	
	public void setMerged(String mergedSeq, String mergedQual, int adjustedAlignmentStart, int adjustedAlignmentEnd) {
		this.mergedSeq = mergedSeq;
		this.mergedQual = mergedQual;		
	}

	public void setMergedSeqAndQual(String mergedSeq, String mergedQual) {
		this.mergedSeq = mergedSeq;
		this.mergedQual = mergedQual;
		
		// Result qual sum
		this.bqSum = -1;
	}
	
	public boolean hasMergedSeq() {
		return this.mergedSeq != null;
	}
	
	public int baseQualSum() {
		if (bqSum < 0) {
			if (hasMergedSeq()) {
				bqSum = SAMRecordUtils.sumBaseQuals(mergedQual);
			} else {
				bqSum = SAMRecordUtils.sumBaseQuals(samRecord);
			}
		}
		
		return bqSum;
	}
	
	public List<Span> getSpanningRegions() {
		
		List<Span> spans = new ArrayList<Span>();
		
		int start = getAdjustedAlignmentStart();
		
		if (samRecord.getReadUnmappedFlag()) {
			spans.add(new Span(start, getAdjustedAlignmentEnd()));
		} else {
			int end = start;
			for (CigarElement elem : samRecord.getCigar().getCigarElements()) {
				switch (elem.getOperator()) {
				case M:
				case S:
				case D:
					end += elem.getLength();
					break;
				case I:
					break;
				case H:
					break;
				case N:
					spans.add(new Span(start, end));
					start = end + elem.getLength();
					end = start;
					break;
				default:
					throw new UnsupportedOperationException("Unhandled cigar operator: " + elem.getOperator() + " in: " + 
							samRecord.getReadName() + " : " + samRecord.getCigarString());
				}
			}
			
			spans.add(new Span(start, end));
		}
		
		return spans;
	}
	
	static class Span {
		int start;
		int end;
		
		public Span(int start, int end) {
			this.start = start;
			this.end = end;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + end;
			result = prime * result + start;
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
			Span other = (Span) obj;
			if (end != other.end)
				return false;
			if (start != other.start)
				return false;
			return true;
		}
	}
}
