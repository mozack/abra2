package abra;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
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
public class ReferenceEvaluator {

	private int readLength;
	private String reference;
	private int[] kmers;
	private String in;
	private String out;
	private String qualities;
	private BufferedWriter include;
	private BufferedWriter exclude;
	private Set<Feature> includeRegions;
	private Set<Feature> excludeRegions;
	private ThreadManager threadManager;
	
	private static final int MAX_NODES = 9000;
	
	public ReferenceEvaluator(int readLength, String reference, int[] kmers, String in, String out, int numThreads) {
		this.readLength = readLength;
		this.reference = reference;
		this.kmers = kmers;
		this.in = in;
		this.out = out;
		
		this.qualities = new String();
		for (int i=0; i<readLength; i++) {
			qualities += "H";
		}
		
		this.threadManager = new ThreadManager(numThreads);
	}
	
	public void run() throws IOException, InterruptedException {
		new NativeLibraryLoader().load(".", NativeLibraryLoader.ABRA, false);
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
		
		System.err.println("Done.");
	}
	
	private void evalRegion(Feature region, String regionBases) {
		boolean shouldInclude = false;
		NativeAssembler assembler = new NativeAssembler();
		StringBuffer readBuf = new StringBuffer((ReAligner.MAX_REGION_LENGTH + 2*readLength) * readLength);
		for (int j=0; j<=regionBases.length() - readLength; j++) {
			readBuf.append("0");  // forward strand only
			String read = regionBases.substring(j, j+readLength);						
			readBuf.append(read);
			readBuf.append(qualities);
		}
		
		String contig = assembler.nativeAssemble(readBuf.toString(), region.getDescriptor(), "eval", 0, 1, (ReAligner.MAX_REGION_LENGTH + 2*readLength)*2, readLength, kmers, 1, 0, .01, 1, MAX_NODES);
		int basesIdx = contig.indexOf('\n') + 1;
		if (basesIdx < contig.length()) {
			String contigBases = contig.substring(basesIdx, contig.length()-1);
			if (regionBases.equals(contigBases)) {
				shouldInclude = true;
			}
		}

		if (shouldInclude) {
			includeRegions.add(region);
		} else {
			excludeRegions.add(region);
		}
	}
	
	private void outputRegions(BufferedWriter writer, Collection<Feature> regions) throws IOException {
		List<Feature> mergedRegions = RegionLoader.collapseRegions(regions, 0);
		
		for (Feature region : mergedRegions) {
			writer.write(region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd() + "\n");
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
		private ReferenceEvaluator evaluator;
		private Feature region;
		private String regionBases;

		public EvalRunnable(ThreadManager threadManager, ReferenceEvaluator evaluator, Feature region, String regionBases) {
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
	
	public static void main(String[] args) throws Exception {
		/*
		NativeLibraryLoader l = new NativeLibraryLoader();
		l.load("/home/lmose/code/abra/target");
		
		int readLength = 100;
		//String reference = "/home/lmose/reference/chr20/20.fa";
		String reference = "/home/lmose/dev/abra/dream/test.fa";
		int[] kmers = new int[] { 43,53,63,73,83 };
		String in = "/home/lmose/dev/abra/dream/include.bed";
		String out = "/home/lmose/dev/abra/dream/exclude.bed";
		*/
		
		if (args.length != 6) {
			System.out.println("ReferenceEvaluator <readLength> <reference> <kmers> <include_bed> <exclude_bed> <num_threads>");
		}
		int readLength = Integer.parseInt(args[0]);
		String reference = args[1];
		int kmers[] = getKmers(args[2]);
		String includeBed = args[3];
		String excludeBed = args[4];
		int threads = Integer.parseInt(args[5]);
		ReferenceEvaluator re = new ReferenceEvaluator(readLength, reference, kmers, includeBed, excludeBed, threads);
		re.run();
	}
}
