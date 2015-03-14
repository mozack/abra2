/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Handles alignment for contigs and short reads.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Aligner {
	
	private String reference;
	private int numThreads;
	private static final int MAX_SMALL_REFERENCE_LINES = 1000000;
	
	public Aligner(String reference, int numThreads) {
		this.reference = reference;
		this.numThreads = numThreads;
	}
	
	public void align(String input, String outputSam, boolean isGapExtensionFavored) throws IOException, InterruptedException {
//		String cmd = "bwa bwasw -t " + numThreads + " -f " + outputSam + " " + reference + " " + input;
		
		String cmd;
		if (!isGapExtensionFavored) {
			cmd = "bwa mem -t " + numThreads + " " + reference + " " + input + " > " + outputSam;
		} else {
			cmd = "bwa mem -A 2 -B 8 -O 12 -L 10 -U 34 -t " + numThreads + " " + reference + " " + input + " > " + outputSam;
		}
		
		runCommand(cmd);
	}

	private void runCommand(String cmd) throws IOException, InterruptedException {
		runCommand(cmd, null, null);
	}
	
	private void runCommand(String cmd, AdjustReadsQueueRunnable stdoutConsumer, ThreadManager threadManager) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		System.out.println("Running: [" + cmd + "]");
		
		long s = System.currentTimeMillis();

		String[] cmds = {
				"bash",
				"-c",
				cmd
			};
		Process proc = Runtime.getRuntime().exec(cmds);
		
//		Thread stdout = new Thread(new CommandOutputConsumer(proc, proc.getInputStream()));
		
		
		Thread stdout = null;
		Thread adjustThread = null;
		ReadInputStreamRunnable isRunnable = null;
		Thread draino = null;
		
		if (stdoutConsumer != null) {
			
			PipedInputStream pis = new PipedInputStream(1000000000);
			PipedOutputStream pos = new PipedOutputStream();
			pos.connect(pis);
			
			// Drain stdout and write to the piped output stream
			draino = new Thread(new Draino(proc.getInputStream(), pos));
			draino.start();
			
			// Read piped input stream and update read queue
			Queue<SAMRecord> queue = new ConcurrentQueue<SAMRecord>();
			stdout = new Thread(new ReadInputStreamRunnable(threadManager, pis, queue));
		 
			// Process read queue content
			stdoutConsumer.setReadQueue(queue);
			adjustThread = new Thread(stdoutConsumer);
			adjustThread.start();
			
		} else {
			stdout = new Thread(new CommandOutputConsumer(proc, proc.getInputStream()));
		}
		
		stdout.start();
		
		Thread stderr = new Thread(new CommandOutputConsumer(proc, proc.getErrorStream()));
		stderr.start();
		
		int ret = proc.waitFor();
		
		stdout.join();
		if (stdoutConsumer != null) {
			stdoutConsumer.setDone();
			adjustThread.join();
		}
		stderr.join();
		
		long e = System.currentTimeMillis();
		
		System.out.println("BWA time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("BWA exited with non-zero return code : [" + ret + "] for command: [" + cmd + "]");
		}
	}
	
	public void shortAlign(String input, String outputSam, AdjustReadsQueueRunnable adjustReadsRunnable, ThreadManager threadManager) throws IOException, InterruptedException {		
		
		//TODO: Just consume bwa output directly?  May allow longer read names.
//		String convert = "bwa samse " + reference + " " + sai + " " + input + " -n 1000 " +
//				"| samtools view -bS -F 0x04 -o " + outputSam + " -";
		
//		String convert = "bwa samse " + reference + " " + sai + " " + input + " -n 1000 " +
//				"| samtools view -bS -o " + outputSam + " -";

		
//		String map = "samtools view " + input + " |  awk '{print \"@\" $1 \"\\n\" $10 \"\\n+\\n\" $11}' | bwa mem -t " + numThreads + " " + reference + " - > " + outputSam;
//		runCommand(map);

		/*
		String sai = outputSam + ".sai";
		
		String aln = "bwa aln " + reference + " " + input + " -f " + sai + " -b -t " + numThreads + " -o 0";
		
		runCommand(aln);
		
		String convert = "bwa samse " + reference + " " + sai + " " + input + " -n 1000 > " + outputSam;
		
		runCommand(convert);
		*/
		
		String map = "bwa aln " + reference + " " + input + " -b -t " + numThreads + " -o 0 | bwa samse " + reference + " - " + input + " -n 1000";
		
		// Redirect stdout to file if no stdout consumer provided.
		if (adjustReadsRunnable == null) {
			map += " > " + outputSam;
		}
		
		runCommand(map, adjustReadsRunnable, threadManager);
	}
	
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
	
	static class CommandOutputConsumer implements Runnable {
		
		private Process proc;
		private InputStream stream;
		
		CommandOutputConsumer(Process proc, InputStream stream) {
			this.proc = proc;
			this.stream = stream;
		}

		@Override
		public void run() {
            InputStreamReader isr = new InputStreamReader(stream);
            BufferedReader br = new BufferedReader(isr);
            String line = null;
            
            try {
	            while ( (line = br.readLine()) != null) {
	            	System.out.println(line);
	            }
	            
	            br.close();
	            isr.close();
            } catch (IOException e) {
            	e.printStackTrace();
            	throw new RuntimeException(e);
            }
            
            System.out.println("Stream thread done.");
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		Aligner a = new Aligner("/datastore/nextgenout2/share/labs/UNCseq/lmose2/mapzilla/bwamem/ref/hg19.fa", 8);
		
		a.align(args[0], args[1], true);
	}
}
