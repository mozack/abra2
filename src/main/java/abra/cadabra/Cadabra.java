package abra.cadabra;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import abra.CompareToReference2;
import abra.Feature;
import abra.Logger;
import abra.ThreadManager;

import abra.cadabra.CadabraProcessor.SampleCall;
import abra.cadabra.CadabraProcessor.SomaticCall;

public class Cadabra {

	private CompareToReference2 c2r;
	
	private Map<String, List<SampleCall>> chromosomeCalls = new HashMap<String, List<SampleCall>>();
	private Map<String, List<SomaticCall>> chromosomeSomaticCalls = new HashMap<String, List<SomaticCall>>();
	
	public void call(String reference, String normalBam, String tumorBam, int numThreads) throws IOException, InterruptedException {
		c2r = new CompareToReference2();
		c2r.init(reference);
		
		outputHeader();
		
		ThreadManager threadManager = new ThreadManager(numThreads);
		
		for (String chromosome : c2r.getChromosomes()) {
			Feature region = new Feature(chromosome, 1, c2r.getChromosomeLength(chromosome));
			CadabraRunnable thread;
			if (normalBam != null) {
				thread = new CadabraRunnable(threadManager, this, normalBam, 
						tumorBam, c2r, region);
			} else {
				thread = new CadabraRunnable(threadManager, this,
						tumorBam, c2r, region);				
			}
			
			threadManager.spawnThread(thread);
		}
		
		threadManager.waitForAllThreadsToComplete();
		
		// Output calls.
		if (normalBam == null) {
			// Simple calling
			for (String chromosome : c2r.getChromosomes()) {
				for (SampleCall call : chromosomeCalls.get(chromosome)) {
					System.out.println(call);
				}
			}
		} else {
			// Somatic calling
			for (String chromosome : c2r.getChromosomes()) {
				for (SomaticCall call : chromosomeSomaticCalls.get(chromosome)) {
					System.out.println(call);
				}
			}
		}
		
		Logger.info("Cadabra done.");
	}
	
	void addCalls(String chromosome, List<SampleCall> calls) {
		Logger.info("Choromosome: %s done.", chromosome);
		synchronized(chromosomeCalls) {
			chromosomeCalls.put(chromosome, calls);
		}
	}
	
	void addSomaticCalls(String chromosome, List<SomaticCall> calls) {
		Logger.info("Choromosome: %s done.", chromosome);
		synchronized(chromosomeSomaticCalls) {
			chromosomeSomaticCalls.put(chromosome, calls);
		}
	}
	
	private void outputHeader() {
		System.out.println("##fileformat=VCFv4.1");
		System.out.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE");
	}
		
	public static void main(String[] args) throws Exception {
//		String normal = "/home/lmose/dev/abra/cadabra/normal_test2.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/tumor_test2.bam";
		
//		String reference = "/home/lmose/reference/chr1/1.fa";
//		String normal = "/home/lmose/dev/abra/cadabra/normal1.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/tumor1.bam";

		
//		String normal = "/home/lmose/dev/abra/cadabra/normal.abra4.sort.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/tumor.abra4.sort.bam";

//		String reference = "/home/lmose/reference/chr1/chr1.fa";
//		String normal = "/home/lmose/dev/abra/cadabra/t2/ntest.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/t2/ttest.bam";

		
//		String reference = "/home/lmose/reference/chr1/chr1.fa";
//		String normal = "/home/lmose/dev/abra/cadabra/ins/ntest.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/ins/ttest.bam";
		
		if (args.length < 3) {
			System.out.println("Usage: java -cp abra.jar abra.cadabra.Cadabra <reference> <bam> <num_threads>");
			System.exit(-1);
		}
		
		String reference = args[0];
		int threads = Integer.parseInt(args[1]);
		String tumor = args[2];
		String normal = null;
		if (args.length > 3) {
			normal = args[3];
		}
		
		new Cadabra().call(reference, normal, tumor, threads);
	}
}
