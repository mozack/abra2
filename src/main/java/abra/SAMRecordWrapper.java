package abra;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SAMRecordWrapper {

	private SAMRecord samRecord;
	private boolean shouldAssemble;
	private boolean shouldFilter;
	private int sampleIdx;
	
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
}
