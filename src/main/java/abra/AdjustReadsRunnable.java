package abra;

import java.io.IOException;

import net.sf.samtools.SAMFileWriter;

public class AdjustReadsRunnable implements Runnable {
	
	private ReAligner realigner;
	private String sortedOriginalReads;
	private String sortedAlignedToContig;
	private SAMFileWriter outputSam;
	private boolean isTightAlignment;
	private CompareToReference2 c2r;
	private String tempDir;
	
//	public AdjustReadsRunnable(ReAligner realigner, String sortedOriginalReads, String sortedAlignedToContig, String outputSam,
//			boolean isTightAlignment) {
	public AdjustReadsRunnable(ReAligner realigner, String sortedAlignedToContig, SAMFileWriter outputSam,
			boolean isTightAlignment, CompareToReference2 c2r, String tempDir) {

		this.realigner = realigner;
		this.sortedOriginalReads = sortedOriginalReads;
		this.sortedAlignedToContig = sortedAlignedToContig;
		this.outputSam = outputSam;
		this.isTightAlignment = isTightAlignment;
		this.c2r = c2r;
		this.tempDir = tempDir;
	}

	@Override
	public void run() {
		try {
			realigner.adjustReads(sortedAlignedToContig, outputSam, isTightAlignment, c2r, tempDir);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
}
