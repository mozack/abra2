package abra;

import java.io.IOException;

public class AdjustReadsRunnable implements Runnable {
	
	private ReAligner realigner;
	private String sortedOriginalReads;
	private String sortedAlignedToContig;
	private String outputSam;
	
	public AdjustReadsRunnable(ReAligner realigner, String sortedOriginalReads, String sortedAlignedToContig, String outputSam) {
		this.realigner = realigner;
		this.sortedOriginalReads = sortedOriginalReads;
		this.sortedAlignedToContig = sortedAlignedToContig;
		this.outputSam = outputSam;
	}

	@Override
	public void run() {
		try {
			realigner.adjustReads(sortedOriginalReads, sortedAlignedToContig, outputSam);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
}
