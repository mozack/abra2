package abra.rna;

import java.util.List;

import net.sf.samtools.SAMRecord;
import abra.AbraRunnable;
import abra.ThreadManager;

public class RnaRegionHandler extends AbraRunnable {
	
	private RnaPoc poc;
	private List<SAMRecord> reads;
	private String id;

	public RnaRegionHandler(ThreadManager threadManager, RnaPoc poc, List<SAMRecord> reads) {
		super(threadManager);
		this.poc = poc;
		this.reads = reads;
		SAMRecord first = reads.get(0);
		SAMRecord last = reads.get(reads.size()-1);
		// NOTE: No guarantee that the last read's alignment end is really the end of the region.
		// Calculating this requires iterating over entire region.
		id = "thread: " + first.getReferenceName() + "_" + first.getAlignmentStart() + "_" + last.getAlignmentEnd() + " reads: " + reads.size();
	}
	
	@Override
	public void go() throws Exception {
		System.out.println("Thread spawned for " + id);
		poc.processReads(reads);
	}
}
