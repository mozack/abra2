package abra;

import java.io.IOException;

public interface StdoutHandler {

	// Buffer up to 1 MB.  Improves performance dramatically when
	// consuming BWA output.
	public static final int MAX_BYTES_TO_BUFFER = 1000000;
	
	/**
	 * Handle processing of stdout from a process.
	 */
	public void process(Process proc) throws IOException;
	
	/**
	 * Handle post-processing of stdout from a process (i.e. thread shutdown / cleanup)
	 */
	public void postProcess() throws InterruptedException;
}
