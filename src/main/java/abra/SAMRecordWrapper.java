package abra;

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
}
