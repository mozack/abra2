package abra.rna;

import java.util.List;

import net.sf.samtools.SAMRecord;
import abra.AbraRunnable;
import abra.ThreadManager;

public class RnaRegionHandler extends AbraRunnable {
	
	private RnaPoc poc;
	private List<SAMRecord> reads;

	public RnaRegionHandler(ThreadManager threadManager, RnaPoc poc, List<SAMRecord> reads) {
		super(threadManager);
		this.poc = poc;
		this.reads = reads;
	}
	
	@Override
	public void go() throws Exception {
		poc.processReads(reads);
	}

}
