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

	public int getAdjustedAlignmentStart() {
		
		//TODO: Adjust unmapped reads by fragment length??
		
		int start = samRecord.getAlignmentStart();
		
		if (samRecord.getCigar().numCigarElements() > 0) {
			CigarElement elem = samRecord.getCigar().getCigarElement(0);
			if (elem.getOperator() == CigarOperator.S) {
				start -= elem.getLength();
				if (start < 1) {
					start = 1;
				}
			}
		}
		
		return start;
	}
	
	public int getAdjustedAlignmentEnd() {
		int end = -1;
		if (samRecord.getReadUnmappedFlag()) {
			// TODO: Pad by fragment length here?
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

		return end;
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
	}
}
