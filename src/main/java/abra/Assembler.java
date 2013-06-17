package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

/**
 * Assembles long contigs from a (relatively small) SAM file.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public interface Assembler {
			
//	public boolean assembleContigs(String inputSam, String output, String prefix, boolean checkForDupes) throws FileNotFoundException, IOException, InterruptedException;
	
	public boolean assembleContigs(List<String> inputFiles, String output, String tempDir, Feature region, String prefix, boolean checkForDupes, ReAligner realigner);
		
//	public void setKmerSize(int kmerSize);
	
	public void setTruncateOutputOnRepeat(boolean truncateOutputOnRepeat);
	
//	public void setMinContigLength(int minContigLength)

//	public void setMinNodeFrequncy(int minNodeFrequncy);
	
	public void setMaxContigs(int maxContigs);
	
//	public void setMinContigRatio(double minContigRatio);
	
	public void setMaxPathsFromRoot(int maxPathsFromRoot);
	
	/*
	static class DepthExceededException extends RuntimeException {

		private int depth;
		
		public DepthExceededException(int depth) {
			this.depth = depth;
		}
		
		public int getDepth() {
			return depth;
		}
	}
	
	static class TooManyPathsFromRootException extends RuntimeException {
		
	}
	
	static class TooManyPotentialContigsException extends RuntimeException {
		
	}
	*/
}
