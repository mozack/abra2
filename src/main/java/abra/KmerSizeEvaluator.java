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
 * Produces BED files indicating genomic windows that lend themselves to assembly.
 *  
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class KmerSizeEvaluator {
	
	public static final int MIN_KMER = 5;

	private int readLength;
	private String reference;
	private String outputFile;
	private String qualities;
	private BufferedWriter output;
	private Set<Feature> outputRegions;
	private ThreadManager threadManager;
	private String regionsBed;
	
	public KmerSizeEvaluator(int readLength, String reference, String outputFile, int numThreads, String regionsBed) {
		this.readLength = readLength;
		this.reference = reference;
		this.outputFile = outputFile;
		this.regionsBed = regionsBed;
		
		this.qualities = new String();
		for (int i=0; i<readLength; i++) {
			qualities += "H";
		}
		
		this.threadManager = new ThreadManager(numThreads);
	}
	
	/*
	public void run() throws IOException, InterruptedException {
//		new NativeLibraryLoader().load(".");
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init8bit(reference);
		
		include = new BufferedWriter(new FileWriter(in, false));
		exclude = new BufferedWriter(new FileWriter(out, false));
		
		for (String chr : c2r.getChromosomes()) {
			includeRegions = Collections.synchronizedSortedSet(new TreeSet<Feature>(new RegionComparator()));
			excludeRegions = Collections.synchronizedSortedSet(new TreeSet<Feature>(new RegionComparator()));
			int i = 0;
			int chromosomeLength = c2r.getReferenceLength(chr);
			while (i < chromosomeLength - ReAligner.MAX_REGION_LENGTH) {
				int regionStart = i;
				int regionStop = i + ReAligner.MAX_REGION_LENGTH;
				int start = Math.max(regionStart - readLength, 0);
				int stop = Math.min(regionStop + readLength, chromosomeLength-1); 
				String regionBases = c2r.getSequence(chr, start+1, stop-start);
				Feature region = new Feature(chr, regionStart, regionStop);
				
				//TODO: Handle other ambiguous bases
				if (!regionBases.contains("N")) {
					threadManager.spawnThread(new EvalRunnable(threadManager, this, region, regionBases));
				} else {
					
					excludeRegions.add(region);
				}
				
				i += ReAligner.REGION_OVERLAP;
			}
			
			threadManager.waitForAllThreadsToComplete();
			
			//TODO: Because assembly regions are overlapped, there is overlap between final include/exclude output
			outputRegions(include, includeRegions);
			outputRegions(exclude, excludeRegions);
		}
		
		include.close();
		exclude.close();
		
		System.out.println("Done.");
	}
	*/
	
	public void run2() throws IOException, InterruptedException {
				
//		new NativeLibraryLoader().load(".");
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init8bit(reference);
		
		List<Feature> regions = ReAligner.getRegions(regionsBed, readLength);
		
		output = new BufferedWriter(new FileWriter(outputFile, false));
		
		outputRegions = Collections.synchronizedSortedSet(new TreeSet<Feature>(new RegionComparator()));
		
		String lastChromosome = "";
		for (Feature region : regions) {
			
			if (!region.getSeqname().equals(lastChromosome)) {
				// Write all entries for this chromosome before moving on to the next
				threadManager.waitForAllThreadsToComplete();
				outputRegions(output, outputRegions);
				outputRegions.clear();
			}
			
			String regionBases = c2r.getSequence(region.getSeqname(), (int) region.getStart()+1-(readLength-1), (int) region.getLength() + (readLength-1));
			
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
		
		System.out.println("Done.");
	}
	
	private void excludeRegion(Feature region) {
		region.setAdditionalInfo("" + (readLength+1) + "\t" + (readLength+1) + "\tN");
	}
	
	private void evalRegion(Feature region, String regionBases) {
		/*
		boolean shouldInclude = false;
		NativeAssembler assembler = new NativeAssembler();
		StringBuffer readBuf = new StringBuffer((ReAligner.MAX_REGION_LENGTH + 2*readLength) * readLength);
		for (int j=0; j<=regionBases.length() - readLength; j++) {
			readBuf.append("0");  // forward strand only
			String read = regionBases.substring(j, j+readLength);						
			readBuf.append(read);
			readBuf.append(qualities);
		}

		int kmer = MIN_KMER;
		while (!shouldInclude && kmer < readLength) {
		
			int kmers[] = new int[] { kmer };
			String contig = assembler.nativeAssemble(readBuf.toString(), region.getDescriptor(), "eval", 0, 1, (ReAligner.MAX_REGION_LENGTH + 2*readLength)*2, readLength, kmers, 1, 0);
			int basesIdx = contig.indexOf('\n') + 1;
			if (basesIdx < contig.length()) {
				String contigBases = contig.substring(basesIdx, contig.length()-1);
				if (regionBases.equals(contigBases)) {
					shouldInclude = true;
				}
			}
			
			if (!shouldInclude) {
				kmer += 2;
			}
		}
		*/
		
		boolean isEditDistanceOK = false;
		int distKmer = MIN_KMER;
		while (!isEditDistanceOK && distKmer < readLength) {
			isEditDistanceOK = isHammingDistanceAtLeast2(regionBases, distKmer);
			if (!isEditDistanceOK) {
				distKmer += 2;
			}
		}

//		region.setAdditionalInfo(String.valueOf(kmer) + "\t" + String.valueOf(distKmer) + "\t.");
		region.setAdditionalInfo(".\t" + String.valueOf(distKmer) + "\t.");
		
//		if (shouldInclude) {
//			region.setAdditionalInfo(String.valueOf(kmer) + "\tINCLUDE");
//		} else {
//			excludeRegion(region);
//		}
		
		outputRegions.add(region);
	}
	
	private void outputRegions(BufferedWriter writer, Collection<Feature> regions) throws IOException {
//		List<Feature> mergedRegions = RegionLoader.collapseRegions(regions, 0);
		
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
	
	private void init(String tempDir) throws IOException {
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

	// Is the hamming distance between all bases in this region at least 2
	private boolean isHammingDistanceAtLeast2(String bases, int k) {
		
		List<String> kmers = new ArrayList<String>();
		for (int i=0; i<=bases.length()-k; i++) {
			String kmer = bases.substring(i, i+k);
			
			for (String kmer2 : kmers) {
				if (!isHammingDistanceAtLeast2(kmer, kmer2)) {
					return false;
				}
			}
			
			kmers.add(kmer);
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
		
//		int readLength = 100;
//		String reference = "/home/lmose/reference/chr20/20.fa";
		//String reference = "/home/lmose/dev/abra/dream/test.fa";
		
		
//		String includeBed = "/home/lmose/dev/abra/dream/round2/include.bed";
//		String excludeBed = "/home/lmose/dev/abra/dream/round2/exclude.bed";
//		int threads = 1;
//		String regionsBed = "/home/lmose/dev/abra/dream/round2/20.orig.bed";
		
		if (args.length != 6) {
			System.out.println("KmerSizeEvaluator <readLength> <reference> <output_bed> <num_threads> <input_bed> <temp_dir>");
		}
		
		int readLength = Integer.parseInt(args[0]);
		String reference = args[1];
		String outputBed = args[2];
		int threads = Integer.parseInt(args[3]);
		String regionsBed = args[4];
		String tempDir = args[5];

		double s = System.currentTimeMillis();
		KmerSizeEvaluator re = new KmerSizeEvaluator(readLength, reference, outputBed, threads, regionsBed);
		re.init(tempDir);
		re.run2();
		double e = System.currentTimeMillis();
		
		double elapsed = (e-s)/1000;
		elapsed /= 60;
		System.out.println("Elapsed mins: " + elapsed);
	}
}
