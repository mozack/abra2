package abra.cadabra;

import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;

public class ReadsAtLocus {

	private String chr;
	private int position;
	private List<SAMRecord> reads;
	
	public ReadsAtLocus(String chr, int position, List<SAMRecord> reads) {
		this.chr = chr;
		this.position = position;
		this.reads = reads;
	}

	public String getChromosome() {
		return chr;
	}

	public int getPosition() {
		return position;
	}

	public List<SAMRecord> getReads() {
		return reads;
	}
	
	public String toString() {
		String s = chr + ":" + position;
		
		for (SAMRecord read : reads) {
			s += "," + read.getReadName();
		}
		
		return s;
	}
	
	public int compareLoci(ReadsAtLocus that, SAMSequenceDictionary dict) {
		int compare = dict.getSequenceIndex(this.getChromosome()) - dict.getSequenceIndex(that.getChromosome());
		if (compare == 0) {
			compare = this.getPosition() - that.getPosition();
		}
		
		return compare;
	}
}
