/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.List;

/**
 * Contig assembly interface.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public interface Assembler {
	
	public List<String> assembleContigs(List<String> inputFiles, String output, String tempDir, Feature region, String prefix, boolean checkForDupes, ReAligner realigner);
	
	public void setTruncateOutputOnRepeat(boolean truncateOutputOnRepeat);
	
	public void setMaxContigs(int maxContigs);
	
	public void setMaxPathsFromRoot(int maxPathsFromRoot);
}
