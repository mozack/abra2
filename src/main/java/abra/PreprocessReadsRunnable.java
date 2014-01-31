package abra;

import static abra.Logger.log;

import net.sf.samtools.SAMFileWriter;

/**
 * Thread entry point for read pre-processing
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class PreprocessReadsRunnable extends AbraRunnable {
	
	private String inputSam;
	private String fastq;
	private CompareToReference2 c2r;
	private SAMFileWriter finalOutputSam;
	private ReAligner reAligner;
	
	public PreprocessReadsRunnable(ThreadManager threadManager, ReAligner reAligner, String inputSam, String fastq, CompareToReference2 c2r, SAMFileWriter finalOutputSam) {
		super(threadManager);
		this.inputSam = inputSam;
		this.fastq = fastq;
		this.c2r = c2r;
		this.finalOutputSam = finalOutputSam;
		this.reAligner = reAligner;
	}

	@Override
	public void go() throws Exception {
		log("Preprocessing original reads for alignment: " + inputSam);
		reAligner.sam2Fastq(inputSam, fastq, c2r, finalOutputSam);
		log("Done preprocessing original reads for alignment: " + inputSam);
	}
}
