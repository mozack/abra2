package abra;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Aligner {
	
	private String reference;
	private int numThreads;
	private static final int MAX_SMALL_REFERENCE_LINES = 1000000;
	
	public Aligner(String reference, int numThreads) {
		this.reference = reference;
		this.numThreads = numThreads;
	}
	
	public void align(String input, String outputSam) throws IOException, InterruptedException {
		String cmd = "bwa bwasw -t " + numThreads + " -f " + outputSam + " " + reference + " " + input;
		
		runCommand(cmd);
	}

	private void runCommand(String cmd) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		System.out.println("Running: [" + cmd + "]");
		
		long s = System.currentTimeMillis();

		String[] cmds = {
				"bash",
				"-c",
				cmd
			};
		Process proc = Runtime.getRuntime().exec(cmds);

//		Process proc;
//		if (inBash) {
//		} else {
//			proc = Runtime.getRuntime().exec(cmd);
//		}
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.out.println("BWA time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("BWA exited with non-zero return code : [" + ret + "] for command: [" + cmd + "]");
		}
	}
	
	public void shortAlign(String input, String outputSam) throws IOException, InterruptedException {
		String sai = outputSam + ".sai";
		
		String aln = "bwa aln " + reference + " " + input + " -f " + sai + " -t " + numThreads + " -o 0";
		
		runCommand(aln);
		
		//TODO: Just consume bwa output directly?  May allow longer read names.
//		String convert = "bwa samse " + reference + " " + sai + " " + input + " -n 1000 " +
//				"| samtools view -bS -F 0x04 -o " + outputSam + " -";
		
		String convert = "bwa samse " + reference + " " + sai + " " + input + " -n 1000 " +
				"| samtools view -bS -o " + outputSam + " -";

		
		runCommand(convert);
	}
	
//	public void smallIndex() throws IOException, InterruptedException {
//		runCommand("bwa index " + reference);
//	}
	
	public void index() throws IOException, InterruptedException {
		
		if (shouldUseSmallIndex()) {
			runCommand("bwa index " + reference);	
		} else {
			runCommand("bwa index -a bwtsw " + reference);
		}
	}
	
	private boolean shouldUseSmallIndex() throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(reference));
		
		try {
			int lines = 0;
			String line = reader.readLine();
			while (line != null) {
				lines++;
				line = reader.readLine();
				
				if (lines >= 1000000) {
					return false;
				}
			}
		} finally {
			reader.close();
		}
		
		return true;
	}
}
