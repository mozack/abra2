package abra.cadabra;

import abra.AbraRunnable;
import abra.CompareToReference2;
import abra.Feature;
import abra.ThreadManager;

public class CadabraRunnable extends AbraRunnable {
	
	private CadabraProcessor processor;
	private Feature region;

	public CadabraRunnable(ThreadManager threadManager, Cadabra cadabra, String bam, 
			CompareToReference2 c2r, Feature region) {
		super(threadManager);
		this.processor = new CadabraProcessor(cadabra, bam, c2r);
		this.region = region;
	}
	
	public CadabraRunnable(ThreadManager threadManager, Cadabra cadabra, 
			String normalBam, String tumorBam, CompareToReference2 c2r, Feature region) {
		super(threadManager);
		this.processor = new CadabraProcessor(cadabra, normalBam, tumorBam, c2r);
		this.region = region;
	}

	@Override
	public void go() throws Exception {
		processor.process(region);
	}

}
