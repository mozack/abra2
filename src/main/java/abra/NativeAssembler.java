package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class NativeAssembler implements Assembler {
	
	private boolean truncateOnRepeat;
	private int maxContigs;
	private int maxPathsFromRoot;
	private Set<String> readIds;

	private native int assemble(String input, String output, String prefix, int truncateOnRepeat, int maxContigs, int maxPathsFromRoot);
//	private native void assemble(String input, String output, String prefix);
	
	static{
        System.loadLibrary("Abra");
	}
	
	private String getIdentifier(SAMRecord read) {
		String id = read.getReadName();
		
		if (read.getReadPairedFlag() && read.getSecondOfPairFlag()) {
			id += "_2";
		}
		
		return id;
	}
	
	public boolean assembleContigs(String input, String output, String prefix, boolean checkForDupes) {
		int count = 0;
		
		readIds = new HashSet<String>();
		
		try {
			SAMFileReader reader = new SAMFileReader(new File(input));
			reader.setValidationStringency(ValidationStringency.SILENT);
			
			String readFile = input + ".reads";
			BufferedWriter writer = new BufferedWriter(new FileWriter(readFile, false));
			
			for (SAMRecord read : reader) {
				// Don't allow same read to be counted twice.
				// TODO: Handle paired end
				if ((!checkForDupes) || (!readIds.contains(getIdentifier(read)))) {
					boolean hasAmbiguousBases = read.getReadString().contains("N");
					Integer numBestHits = (Integer) read.getIntegerAttribute("X0");
					boolean hasAmbiguousInitialAlignment = numBestHits != null && numBestHits > 1;
					//TODO: Stampy ambiguous read (mapq < 4)
					
					if (!hasAmbiguousBases && !hasAmbiguousInitialAlignment) {
						if (!checkForDupes) {
							readIds.add(getIdentifier(read));
						}
						writer.write(read.getReadString() + "\n");
					}
				}
			}
			
			readIds = null;
			
			writer.close();
			reader.close();
			
			count = assemble(
					readFile,
					output, 
					prefix, 
					truncateOnRepeat ? 1 : 0,
					maxContigs,
					maxPathsFromRoot);

		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		
		return count != 0;
	}

	public boolean isTruncateOnRepeat() {
		return truncateOnRepeat;
	}

	public void setTruncateOutputOnRepeat(boolean truncateOnRepeat) {
		this.truncateOnRepeat = truncateOnRepeat;
	}

	public int getMaxContigs() {
		return maxContigs;
	}

	public void setMaxContigs(int maxContigs) {
		this.maxContigs = maxContigs;
	}

	public int getMaxPathsFromRoot() {
		return maxPathsFromRoot;
	}

	public void setMaxPathsFromRoot(int maxPathsFromRoot) {
		this.maxPathsFromRoot = maxPathsFromRoot;
	}
	
	public static void run(String input, String output) {
		NativeAssembler assem = new NativeAssembler();
		assem.setTruncateOutputOnRepeat(false);
		assem.setMaxContigs(2000000);
		assem.setMaxPathsFromRoot(5000);
		
		assem.assembleContigs(input, output, "contig", true);
		
	}
	
	public static void main(String[] args) {
		NativeAssembler assem = new NativeAssembler();
		assem.setTruncateOutputOnRepeat(false);
		assem.setMaxContigs(2000000);
		assem.setMaxPathsFromRoot(5000);
		/*
		assem.assembleContigs(args[0], args[1], "contig");
		*/
		
//		for (int i=0; i<10; i++) {
//			run(args[0], args[1] + "_" + i);
//		}
		
		run(args[0], args[1]);
		
//		assem.assembleContigs("/home/lmose/code/abra/src/main/c/1810_reads.txt",
//				"/home/lmose/code/abra/src/main/c/1810.fa", "bar");
	}
}
