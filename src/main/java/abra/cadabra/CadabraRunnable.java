package abra.cadabra;

import abra.AbraRunnable;
import abra.CompareToReference2;
import abra.Feature;
import abra.ThreadManager;

public class CadabraRunnable extends AbraRunnable {
	
	private CadabraProcessor processor;
	private Feature region;

	public CadabraRunnable(ThreadManager threadManager, Cadabra cadabra, CadabraOptions options, 
			CompareToReference2 c2r, Feature region) {
		super(threadManager);
		this.processor = new CadabraProcessor(cadabra, options, c2r);
		this.region = region;
	}

	@Override
	public void go() throws Exception {
		processor.process(region);
	}

}
