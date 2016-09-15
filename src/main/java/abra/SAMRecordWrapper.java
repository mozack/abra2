package abra;

import htsjdk.samtools.SAMRecord;

public class SAMRecordWrapper {

	private SAMRecord samRecord;
	private boolean shouldAssemble;
	
	public SAMRecordWrapper(SAMRecord record, boolean shouldAssemble) {
		this.samRecord = record;
		this.shouldAssemble = shouldAssemble;
	}

	public SAMRecord getSamRecord() {
		return samRecord;
	}

	public boolean shouldAssemble() {
		return shouldAssemble;
	}
}
