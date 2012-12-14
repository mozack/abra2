package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class NativeAssembler implements Assembler {
	
	private boolean truncateOnRepeat;
	private int maxContigs;
	private int maxPathsFromRoot;

	private native int assemble(String input, String output, String prefix, int truncateOnRepeat, int maxContigs, int maxPathsFromRoot);
//	private native void assemble(String input, String output, String prefix);
	
	static{
        System.loadLibrary("Abra");
	}
	
	public boolean assembleContigs(String input, String output, String prefix) {
		int count = 0;
		
		try {
			SAMFileReader reader = new SAMFileReader(new File(input));
			reader.setValidationStringency(ValidationStringency.SILENT);
			
			String readFile = input + ".reads";
			BufferedWriter writer = new BufferedWriter(new FileWriter(readFile, false));
			
			for (SAMRecord read : reader) {
				boolean hasAmbiguousBases = read.getReadString().contains("N");
				Integer numBestHits = (Integer) read.getIntegerAttribute("X0");
				boolean hasAmbiguousInitialAlignment = numBestHits != null && numBestHits > 1;
				
				if (!hasAmbiguousBases && !hasAmbiguousInitialAlignment) {
					writer.write(read.getReadString() + "\n");
				}
			}
			
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
		
		assem.assembleContigs(input, output, "contig");
		
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
