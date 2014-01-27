package abra;

import static abra.Logger.log;

import java.io.IOException;

import net.sf.samtools.SAMFileWriter;

/**
 * Thread entry point for read pre-processing
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class PreprocessReadsRunnable implements Runnable {
	
	private ReAligner realigner;
	private String inputSam;
	private String fastq;
	private CompareToReference2 c2r;
	private SAMFileWriter finalOutputSam;
	
	public PreprocessReadsRunnable(ReAligner realigner, String inputSam, String fastq, CompareToReference2 c2r, SAMFileWriter finalOutputSam) {
		this.realigner = realigner;
		this.inputSam = inputSam;
		this.fastq = fastq;
		this.c2r = c2r;
		this.finalOutputSam = finalOutputSam;
	}

	@Override
	public void run() {
		log("Preprocessing original reads for alignment: " + inputSam);
		try {
			realigner.sam2Fastq(inputSam, fastq, c2r, finalOutputSam);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		} finally {
			realigner.removeThread(this);
		}
		log("Done preprocessing original reads for alignment: " + inputSam);
	}
}
