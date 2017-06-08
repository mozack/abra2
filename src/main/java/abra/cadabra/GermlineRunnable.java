package abra.cadabra;

import abra.AbraRunnable;
import abra.CompareToReference2;
import abra.Feature;
import abra.ThreadManager;

public class GermlineRunnable extends AbraRunnable {
	
	private GermlineProcessor processor;
	private Feature region;

	public GermlineRunnable(ThreadManager threadManager, Germline germline, String bam, 
			CompareToReference2 c2r, Feature region) {
		super(threadManager);
		this.processor = new GermlineProcessor(germline, bam, c2r);
		this.region = region;
	}
	
	public GermlineRunnable(ThreadManager threadManager, Germline germline, 
			String normalBam, String tumorBam, CompareToReference2 c2r, Feature region) {
		super(threadManager);
		this.processor = new GermlineProcessor(germline, normalBam, tumorBam, c2r);
		this.region = region;
	}

	@Override
	public void go() throws Exception {
		processor.process(region);
	}

}
