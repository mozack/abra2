package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/**
 * Produces BED files indicating genomic windows and kmer sizes that lend themselves to assembly.
 *  
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class KmerSizeEvaluator {
	
	public static final int MIN_KMER = 9;

	private int readLength;
	private CompareToReference2 c2r;
	private String outputFile;
	private BufferedWriter output;
	private Set<Feature> outputRegions;
	private ThreadManager threadManager;
	private String regionsBed;
	
	public KmerSizeEvaluator(int readLength, CompareToReference2 c2r, String outputFile, int numThreads, String regionsBed) {
		this.readLength = readLength;
		this.c2r = c2r;
		this.outputFile = outputFile;
		this.regionsBed = regionsBed;
		this.threadManager = new ThreadManager(numThreads);
	}
	
	public KmerSizeEvaluator() {
	}
	
	private String getBases(Feature region, CompareToReference2 c2r) {
		return c2r.getSequence(region.getSeqname(), (int) region.getStart()+1-(readLength-1), (int) region.getLength() + (readLength*2-2));
	}
		
	public void run() throws IOException, InterruptedException {
				
//		new NativeLibraryLoader().load(".");		
		List<Feature> regions = ReAligner.getRegions(regionsBed, readLength, false);
		
		output = new BufferedWriter(new FileWriter(outputFile, false));
		
		outputRegions = Collections.synchronizedSortedSet(new TreeSet<Feature>(new RegionComparator()));
		
		String lastChromosome = "";
		for (Feature region : regions) {
			
			if (!region.getSeqname().equals(lastChromosome)) {
				// Write all entries for this chromosome before moving on to the next
				threadManager.waitForAllThreadsToComplete();
				outputRegions(output, outputRegions);
				outputRegions.clear();
				threadManager = new ThreadManager(threadManager.getNumThreads());
			}
			
			String regionBases = getBases(region, c2r);
			
			//TODO: Handle other ambiguous bases
			if (!regionBases.contains("N")) {
				threadManager.spawnThread(new EvalRunnable(threadManager, this, region, regionBases));
			} else {
				excludeRegion(region);
				outputRegions.add(region);
			}
			
			lastChromosome = region.getSeqname();
		}
			
		threadManager.waitForAllThreadsToComplete();
		
		// Because assembly regions are overlapped, there is overlap between final include/exclude output
		outputRegions(output, outputRegions);
		
		output.close();
		
		System.err.println("Done.");
	}
	
	private void excludeRegion(Feature region) {
		region.setAdditionalInfo((readLength+1) + "\tN");
	}
	
	public int identifyMinKmer(int readLength, CompareToReference2 c2r, List<Feature> regions) {
		List<String> regionBases = new ArrayList<String>();
		
		for (Feature region : regions) {
			if (region.getLength() > MIN_KMER) {
				regionBases.add(getBases(region, c2r));
			}
		}
		
		boolean isEditDistanceOK = false;
		int distKmer = MIN_KMER;
		while (!isEditDistanceOK && distKmer < readLength) {
			isEditDistanceOK = isHammingDistanceAtLeast2(regionBases, distKmer);
			if (!isEditDistanceOK) {
				distKmer += 2;
			}
		}

		return distKmer;
	}
	
	private void evalRegion(Feature region, String regionBases) {	

		List<Feature> regionList = new ArrayList<Feature>();
		regionList.add(region);
		int distKmer = identifyMinKmer(readLength, c2r, regionList);
		
		region.setAdditionalInfo(String.valueOf(distKmer) + "\t.");
				
		outputRegions.add(region);
	}
	
	private void outputRegions(BufferedWriter writer, Collection<Feature> regions) throws IOException {
		
		for (Feature region : regions) {
			String output = region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd();
			if (region.getAdditionalInfo() != null) {
				output += "\t" + region.getAdditionalInfo();
			}
			output += "\n";
			writer.write(output);
		}
	}
	
	// Basic Region comparator.  Only considers start position, so must not
	// be used across chromosomes.
	static class RegionComparator implements Comparator<Feature> {

		@Override
		public int compare(Feature region1, Feature region2) {
			return (int) (region1.getStart() - region2.getStart());
		}
	}
	
	static class EvalRunnable extends AbraRunnable {
		private KmerSizeEvaluator evaluator;
		private Feature region;
		private String regionBases;

		public EvalRunnable(ThreadManager threadManager, KmerSizeEvaluator evaluator, Feature region, String regionBases) {
			super(threadManager);
			this.evaluator = evaluator;
			this.region = region;
			this.regionBases = regionBases;
		}

		@Override
		public void go() throws Exception {
			evaluator.evalRegion(region, regionBases);
		}
	}
	
	static int[] getKmers(String str) {
		String[] strings = str.split(",");
		int[] kmers = new int[strings.length];
		for (int i=0; i<strings.length; i++) {
			kmers[i] = Integer.parseInt(strings[i]);
		}
		return kmers;
	}
	
	// Is the hamming distance between all bases in this region at least 2
	private boolean isHammingDistanceAtLeast2(List<String> basesList, int k) {
		
		List<String> kmers = new ArrayList<String>();
		
		for (String bases : basesList) {
			for (int i=0; i<=bases.length()-k; i++) {
				String kmer = bases.substring(i, i+k);
				
				for (String kmer2 : kmers) {
					if (!isHammingDistanceAtLeast2(kmer, kmer2)) {
						return false;
					}
				}
				
				kmers.add(kmer);
			}
		}
		
		return true;
	}
	
	private boolean isHammingDistanceAtLeast2(String kmer1, String kmer2) {
		int dist = 0;
		for (int i=0; i<kmer1.length(); i++) {
			if (kmer1.charAt(i) != kmer2.charAt(i)) {
				dist += 1;
				
				if (dist >=2) {
					return true;
				}
			}
		}
		
		return false;
	}
	
	public static void main(String[] args) throws Exception {
		
//		NativeLibraryLoader l = new NativeLibraryLoader();
//		l.load("/home/lmose/code/abra/target");
//		
//		int readLength = 100;
//		String reference = "/home/lmose/reference/chr20/20.fa";
//		//String reference = "/home/lmose/dev/abra/dream/test.fa";
//		
//		
//		int threads = 1;
//		String outputBed = "/home/lmose/dev/abra/dream/round2/output.bed";
//		String regionsBed = "/home/lmose/dev/abra/dream/round2/20.orig.bed";
		
		if (args.length != 5) {
			System.err.println("KmerSizeEvaluator <readLength> <reference> <output_bed> <num_threads> <input_bed>");
			System.exit(-1);
		}
		
		int readLength = Integer.parseInt(args[0]);
		String reference = args[1];
		String outputBed = args[2];
		int threads = Integer.parseInt(args[3]);
		String regionsBed = args[4];

		double s = System.currentTimeMillis();
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(reference);
		
		KmerSizeEvaluator re = new KmerSizeEvaluator(readLength, c2r, outputBed, threads, regionsBed);
		re.run();
		double e = System.currentTimeMillis();
		
		double elapsed = (e-s)/1000;
		elapsed /= 60;
		System.err.println("Elapsed mins: " + elapsed);
	}
}
