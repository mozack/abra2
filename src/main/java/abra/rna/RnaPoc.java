package abra.rna;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import abra.NativeAssembler;
import abra.NativeLibraryLoader;
import abra.ThreadManager;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class RnaPoc {
	
	//public static int MAX_READ_GAP = 1000000;
	public static int MAX_READ_GAP = 50000;
	
	private BufferedWriter contigWriter;
	
	private ThreadManager threadManager;
	
	private void init(String tempDir) {
		File workingDir = new File(tempDir);
		if (workingDir.exists()) {
			if (!workingDir.delete()) {
				throw new IllegalStateException("Unable to delete: " + tempDir);
			}
		}

		if (!workingDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir);
		}
						
		new NativeLibraryLoader().load(tempDir);

	}

	public void run(String input, String output, String temp, int numThreads) throws IOException, InterruptedException {
		
		init(temp);
		
		this.threadManager = new ThreadManager(numThreads);
		
		contigWriter = new BufferedWriter(new FileWriter(output, false));
		
		List<SAMRecord> currReads = new ArrayList<SAMRecord>();
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		int prevMaxEnd = -1;
		
		SAMRecord lastRead = null;

		for (SAMRecord read : reader) {
			if (read.getMappingQuality() > 0) {
				if (lastRead == null || !lastRead.getReferenceName().equals(read.getReferenceName()) || (read.getAlignmentStart()-prevMaxEnd) < MAX_READ_GAP) {
					currReads.add(read);
				} else {
//					processReads(currReads);
					spawnProcessingThread(currReads);
					currReads = new ArrayList<SAMRecord>();
					currReads.add(read);
				}
				
				if (read.getAlignmentEnd() > prevMaxEnd || !lastRead.getReferenceName().equals(read.getReferenceName())) {
					prevMaxEnd = read.getAlignmentEnd();
				}
				
				lastRead = read;
			}
		}
		
		threadManager.waitForAllThreadsToComplete();
		
		reader.close();
		contigWriter.close();
	}
	
	private void spawnProcessingThread(List<SAMRecord> reads) {
		RnaRegionHandler handler = new RnaRegionHandler(threadManager, this, reads);
		threadManager.spawnThread(handler);
	}
	
	void processReads(List<SAMRecord> reads) throws IOException {

		NativeAssembler assem = newAssembler();
		
		String contigs = assem.simpleAssemble(reads);
		
		if (!contigs.equals("<ERROR>") && !contigs.equals("<REPEAT>") && !contigs.isEmpty()) {
			appendContigs(contigs);
		}
	}
	
	private synchronized void appendContigs(String contigs) throws IOException {
		contigWriter.write(contigs);
	}
	
	private NativeAssembler newAssembler() {
		NativeAssembler assem = new NativeAssembler();

		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(1000);

		assem.setMaxPathsFromRoot(5000000);
		assem.setReadLength(75);
		assem.setKmer(new int[] { 17, 27, 37, 47 });
		assem.setMinKmerFrequency(2);
		assem.setMinBaseQuality(40);
		
		// The following params not used
		assem.setMinReadCandidateFraction(0);
		assem.setMaxAverageDepth(0);
		assem.setShouldSearchForSv(false);
		assem.setAverageDepthCeiling(0);

		return assem;		
	}
	
	public static void main(String[] args) throws Exception {
		RnaPoc poc = new RnaPoc();
		
		poc.run(args[0], args[1], args[2], Integer.parseInt(args[3]));
	}
}
