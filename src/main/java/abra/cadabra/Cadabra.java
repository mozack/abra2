package abra.cadabra;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import abra.CompareToReference2;
import abra.Feature;
import abra.Logger;
import abra.ThreadManager;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class Cadabra {

	private CompareToReference2 c2r;
	
	private Map<String, List<String>> chromosomeCalls;

	public void callSomatic(String reference, String normal, String tumor, int numThreads) throws IOException, InterruptedException {
		c2r = new CompareToReference2();
		c2r.init(reference);
		
//		if (target != null) {
//			String[] fields = target.split(":|-");
//			String chromosome = fields[0];
//			long startPos = Long.valueOf(fields[1]);
//			long endPos = Long.valueOf(fields[2]);
//			region = new Feature(chromosome, startPos, endPos);
//		}

		outputHeader();
		
		ThreadManager threadManager = new ThreadManager(numThreads);
		
		for (String chromosome : c2r.getChromosomes()) {
			Feature region = new Feature(chromosome, 1, c2r.getChromosomeLength(chromosome));
			CadabraRunnable thread = new CadabraRunnable(threadManager, this, normal, tumor, 
				c2r, region);
			
			threadManager.spawnThread(thread);
		}
		
		threadManager.waitForAllThreadsToComplete();
		
		// Output calls.
		for (String chromosome : c2r.getChromosomes()) {
			for (String call : chromosomeCalls.get(chromosome)) {
				System.out.println(call);
			}
		}
		
		Logger.info("Cadabra done.");
	}
	
	void addCalls(String chromosome, List<String> calls) {
		Logger.info("Choromosome: %s done.", chromosome);
		synchronized(calls) {
			chromosomeCalls.put(chromosome, calls);
		}
	}
	
	private void outputHeader() {
		System.out.println("##fileformat=VCFv4.1");
		System.out.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR");
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
			System.out.println("Usage: java -cp abra.jar abra.cadabra.Cadabra <reference> <normal_bam> <tumor_bam> <num_threads>");
			System.exit(-1);
		}
		
		String reference = args[0];
		String normal = args[1];
		String tumor = args[2];
		int threads = Integer.parseInt(args[3]);
		
		new Cadabra().callSomatic(reference, normal, tumor, threads);
	}
}
