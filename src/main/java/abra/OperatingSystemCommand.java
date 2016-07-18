/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;

/**
 * Handles spawning of programs.  Likely to be replaced with JNI.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class OperatingSystemCommand {

	public static void runCommand(String cmd) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		System.err.println("Running: [" + cmd + "]");
		
		long s = System.currentTimeMillis();
		
		Process proc = Runtime.getRuntime().exec(cmd);
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.err.println("cmd time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("cmd exited with non-zero return code : [" + ret + "] for command: [" + cmd + "]");
		}
		
		proc.destroy();
	}

	public static void runCommand(String[] cmd) throws IOException, InterruptedException {
		
		//String cmd = "bwa bwasw -f " + outputSam + " " + reference + " " + input;
		StringBuffer cmdStr = new StringBuffer();
		for (String substring : cmd) {
			cmdStr.append(substring);
			cmdStr.append(" ");
		}
		
		System.err.println("Running: [" + cmdStr + "]");
		
		long s = System.currentTimeMillis();
		
		Process proc = Runtime.getRuntime().exec(cmd);
		
		//TODO: Catch InterruptedException ?
		//TODO: Capture stderr
		int ret = proc.waitFor();
		
		long e = System.currentTimeMillis();
		
		System.err.println("cmd time: " + (e-s)/1000 + " seconds.");
		
		if (ret != 0) {
			throw new RuntimeException("cmd exited with non-zero return code : [" + ret + "] for command: [" + cmdStr + "]");
		}
		
		proc.destroy();
	}
}
