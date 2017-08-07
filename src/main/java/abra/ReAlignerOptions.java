/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import joptsimple.OptionParser;

/**
 * Manages ABRA command line options
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAlignerOptions extends Options {
	
	private static final String INPUT_SAM = "in";
	private static final String OUTPUT_SAM = "out";
	private static final String REFERENCE = "ref";
	private static final String TARGET_REGIONS = "targets";
	private static final String TARGET_REGIONS_WITH_KMERS = "target-kmers";
	private static final String KMER_SIZE = "kmer";
	private static final String MIN_NODE_FREQUENCY = "mnf";
	private static final String MIN_CONTIG_LENGTH = "mcl";
	private static final String MIN_MAPQ = "mapq";
	private static final String NUM_THREADS = "threads";
	private static final String SINGLE_END = "single";
	private static final String MIN_BASE_QUALITY = "mbq";
	private static final String MIN_READ_CANDIDATE_FRACTION = "rcf";
	private static final String MAX_AVERAGE_REGION_DEPTH = "mad";
	private static final String MIN_EDGE_RATIO = "mer";
	private static final String MAX_NODES = "maxn";
	private static final String SKIP_ASSEMBLY = "sa";
	private static final String SKIP_SOFT_CLIP = "ssc";
	private static final String SOFT_CLIP = "sc";
	private static final String SKIP_OBS_INDELS = "sobs";
	private static final String JUNCTIONS = "junctions";
	private static final String LOG_LEVEL = "log";
	private static final String CONTIG_FILE	 = "contigs";
	private static final String GTF_JUNCTIONS = "gtf";
	private static final String SG_ALIGNER_SCORING = "sga";
	private static final String MAX_CACHED_READS = "mcr";
	private static final String KEEP_TMP = "keep-tmp";
	private static final String TMP_DIR = "tmpdir";
	private static final String CONSENSUS_SEQ = "cons";
	private static final String MAX_MISMATCH_RATE = "mmr";
	private static final String WINDOW_SIZE = "ws";
	private static final String MAX_READS_IN_REGION = "mrr";
	private static final String COMPRESSION_LEVEL = "cl";
	private static final String CONTIG_ANCHOR = "ca";
	private static final String NO_SORT = "nosort";
	private static final String MAX_READ_MOVE_DISTANCE = "dist";
	private static final String CHROMOSOMES_TO_SKIP = "skip";
	private static final String MAX_ASSEMBLED_CONTIGS = "mac";
	private static final String SKIP_UNMAPPED_ASSEMBLY_TRIGGER = "sua";
	private static final String UNSET_DUPLICATES = "undup";
	private static final String INPUT_VCF = "in-vcf";
	private static final String INDEX = "index";
	private static final String NO_GKL = "no-gkl";
	
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required list of input sam or bam file(s) separated by comma").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM, "Required list of output sam or bam file(s) separated by comma").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE, "Genome reference location").withRequiredArg().ofType(String.class);
            parser.accepts(TARGET_REGIONS, "BED file containing target regions").withRequiredArg().ofType(String.class);
            parser.accepts(TARGET_REGIONS_WITH_KMERS, "BED-like file containing target regions with per region kmer sizes in 4th column").withRequiredArg().ofType(String.class);
            parser.accepts(KMER_SIZE, "Optional assembly kmer size(delimit with commas if multiple sizes specified)").withOptionalArg().ofType(String.class);
            parser.accepts(MIN_NODE_FREQUENCY, "Assembly minimum node frequency").withRequiredArg().ofType(Integer.class).defaultsTo(1);
            parser.accepts(MIN_CONTIG_LENGTH, "Assembly minimum contig length").withOptionalArg().ofType(Integer.class).defaultsTo(-1);
            parser.accepts(NUM_THREADS, "Number of threads").withRequiredArg().ofType(Integer.class).defaultsTo(4);
            parser.accepts(MIN_MAPQ, "Minimum mapping quality for a read to be used in assembly and be eligible for realignment").withOptionalArg().ofType(Integer.class).defaultsTo(20);
            parser.accepts(SINGLE_END, "Input is single end");
            parser.accepts(MIN_BASE_QUALITY, "Minimum base quality for inclusion in assembly.  This value is compared against the sum of base qualities per kmer position").withOptionalArg().ofType(Integer.class).defaultsTo(20);
            parser.accepts(MIN_READ_CANDIDATE_FRACTION, "Minimum read candidate fraction for triggering assembly").withRequiredArg().ofType(Double.class).defaultsTo(.01);
            parser.accepts(MAX_AVERAGE_REGION_DEPTH, "Regions with average depth exceeding this value will be downsampled").withRequiredArg().ofType(Integer.class).defaultsTo(1000);
            parser.accepts(MIN_EDGE_RATIO, "Min edge pruning ratio.  Default value is appropriate for relatively sensitive somatic cases.  May be increased for improved speed in germline only cases.").withRequiredArg().ofType(Double.class).defaultsTo(.01);
            parser.accepts(MAX_NODES, "Maximum pre-pruned nodes in regional assembly").withOptionalArg().ofType(Integer.class).defaultsTo(150000);
            parser.accepts(SKIP_ASSEMBLY, "Skip assembly");
            parser.accepts(JUNCTIONS, "Splice junctions definition file").withRequiredArg().ofType(String.class);
            parser.accepts(LOG_LEVEL, "Logging level (trace,debug,info,warn,error)").withRequiredArg().ofType(String.class).defaultsTo("info");
            parser.accepts(CONTIG_FILE, "Optional file to which assembled contigs are written").withRequiredArg().ofType(String.class);
            parser.accepts(GTF_JUNCTIONS, "GTF file defining exons and transcripts").withRequiredArg().ofType(String.class);
            parser.accepts(SKIP_SOFT_CLIP, "Skip usage of soft clipped sequences as putative contigs");
            parser.accepts(SOFT_CLIP, "Soft clip contig args [max_contigs,min_base_qual,frac_high_qual_bases,min_soft_clip_len]").withRequiredArg().ofType(String.class).defaultsTo("16,13,80,15");
            parser.accepts(SG_ALIGNER_SCORING, "Scoring used for contig alignments (match, mismatch_penalty, gap_open_penalty, gap_extend_penalty)").withRequiredArg().ofType(String.class).defaultsTo("8,32,48,1");
            parser.accepts(MAX_CACHED_READS, "Max number of cached reads per sample per thread").withRequiredArg().ofType(Integer.class).defaultsTo(1000000);
            parser.accepts(SKIP_OBS_INDELS, "Do not use observed indels in original alignments to generate contigs");
            parser.accepts(KEEP_TMP, "Do not delete the temporary directory");
            parser.accepts(TMP_DIR, "Set the temp directory (overrides java.io.tmpdir)").withRequiredArg().ofType(String.class);
            parser.accepts(CONSENSUS_SEQ, "Use positional consensus sequence when aligning high quality soft clipping");
            parser.accepts(MAX_MISMATCH_RATE, "Max allowed mismatch rate when mapping reads back to contigs").withRequiredArg().ofType(Double.class).defaultsTo(.05);
            parser.accepts(WINDOW_SIZE, "Processing window size and overlap (size,overlap)").withRequiredArg().ofType(String.class).defaultsTo("400,200");
            parser.accepts(MAX_READS_IN_REGION, "Regions containing more reads than this value are not processed.  Use -1 to disable.").withRequiredArg().ofType(Integer.class).defaultsTo(1000000);
            parser.accepts(COMPRESSION_LEVEL, "Compression level of output bam file(s)").withRequiredArg().ofType(Integer.class).defaultsTo(5);
            parser.accepts(CONTIG_ANCHOR, "Contig anchor [M_bases_at_contig_edge, max_mismatches_near_edge]").withRequiredArg().ofType(String.class).defaultsTo("10,2");
            parser.accepts(NO_SORT, "Do not attempt to sort final output");
            parser.accepts(MAX_READ_MOVE_DISTANCE, "Max read move distance").withRequiredArg().ofType(Integer.class).defaultsTo(1000);
            parser.accepts(CHROMOSOMES_TO_SKIP, "If no target specified, skip realignment of chromosomes matching specified regex.  Skipped reads are output without modification.  Specify none to disable.").withRequiredArg().ofType(String.class).defaultsTo(ChromosomeRegex.DEFAULT_SKIP_REGEX);
            parser.accepts(MAX_ASSEMBLED_CONTIGS, "Max assembled contigs").withRequiredArg().ofType(Integer.class).defaultsTo(64);
            parser.accepts(SKIP_UNMAPPED_ASSEMBLY_TRIGGER, "Do not use unmapped reads anchored by mate to trigger assembly.  These reads are still eligible to contribute to assembly");
            parser.accepts(UNSET_DUPLICATES, "Unset duplicate flag");
            parser.accepts(INPUT_VCF, "VCF containing known (or suspected) variant sites.  Very large files should be avoided.").withRequiredArg().ofType(String.class);
            parser.accepts(INDEX, "Enable BAM index generation when outputting sorted alignments (may require additonal memory)");
            parser.accepts(NO_GKL, "Disable GKL Intel Deflater");
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
		isValid = true;
		
		if (!getOptions().hasArgument(INPUT_SAM)) {
			isValid = false;
			System.err.println("Missing required input SAM/BAM file");
		}

		if (!getOptions().hasArgument(OUTPUT_SAM)) {
			isValid = false;
			System.err.println("Missing required input SAM/BAM file");
		}
		
		if (getInputFiles().length != getOutputFiles().length) {
			System.err.println("Number of input files must equal number of output files");
		}
		
		if (!getOptions().hasArgument(REFERENCE)) {
			isValid = false;
			System.err.println("Missing required reference");
		}
		
		if (getOptions().hasArgument(TARGET_REGIONS) && getOptions().hasArgument(TARGET_REGIONS_WITH_KMERS)) {
			isValid = false;
			System.err.println("Please specifiy only one of: " + TARGET_REGIONS + ", " + TARGET_REGIONS_WITH_KMERS);
		}		
				
		if ((getOptions().hasArgument(NUM_THREADS) && (Integer) getOptions().valueOf(NUM_THREADS) < 1)) {
			isValid = false;
			System.err.println("Num threads must be greater than zero.");
		}
		
        if (!isValid) {
            printHelp();
        }
	}
	
	public String[] getInputFiles() {
		String[] files = new String[0];
		String sams = (String) getOptions().valueOf(INPUT_SAM);
		if (sams != null) {
			files = sams.split(",");
		}
		return files;
	}
	
	public String[] getOutputFiles() {
		String[] files = new String[0];
		String sams = (String) getOptions().valueOf(OUTPUT_SAM);
		if (sams != null) {
			files = sams.split(",");
		}
		return files;
	}
	
	public String getChromosomesToSkipRegex() {
		return (String) getOptions().valueOf(CHROMOSOMES_TO_SKIP);
	}
		
	public String getReference() {
		return (String) getOptions().valueOf(REFERENCE);
	}
	
	public String getTargetRegionFile() {
		
		String file = null;
		
		if (getOptions().hasArgument(TARGET_REGIONS_WITH_KMERS)) {
			file = (String) getOptions().valueOf(TARGET_REGIONS_WITH_KMERS);
		} else {
			file = (String) getOptions().valueOf(TARGET_REGIONS);
		}
		
		return file;
	}
	
	public String getJunctionFile() {
		String file = null;
		
		if (getOptions().hasArgument(JUNCTIONS)) {
			file = (String) getOptions().valueOf(JUNCTIONS);
		}
		
		return file;		
	}
	
	public String getGtfJunctionFile() {
		String file = null;
		
		if (getOptions().hasArgument(GTF_JUNCTIONS)) {
			file = (String) getOptions().valueOf(GTF_JUNCTIONS);
		}
		
		return file;		
	}

	public boolean hasPresetKmers() {
		return getOptions().hasArgument(TARGET_REGIONS_WITH_KMERS);
	}
	
	public int[] getKmerSizes() {
		int[] kmers;
		if (getOptions().has(KMER_SIZE)) {
			String[] kmerStr = ((String) getOptions().valueOf(KMER_SIZE)).split(",");
			
			kmers = new int[kmerStr.length];
			for (int i=0; i<kmerStr.length; i++) {
				kmers[i] = Integer.parseInt(kmerStr[i]);
			}
		} else {
			kmers = new int[0];
		}
		
		return kmers;
	}
		
	public int getMinNodeFrequency() {
		return (Integer) getOptions().valueOf(MIN_NODE_FREQUENCY);
	}
	
	public int getMinContigLength() {
		return (Integer) getOptions().valueOf(MIN_CONTIG_LENGTH);
	}
	
	public int getMaxRealignDist() {
		return (Integer) getOptions().valueOf(MAX_READ_MOVE_DISTANCE);
	}
	
	public int getMaxAssembledContigs() {
		return (Integer) getOptions().valueOf(MAX_ASSEMBLED_CONTIGS);
	}
	
	public int getNumThreads() {
		return getOptions().hasArgument(NUM_THREADS) ? (Integer) getOptions().valueOf(NUM_THREADS) : 4;
	}
	
	public int getMaxReadsInRegion() {
		return (Integer) getOptions().valueOf(MAX_READS_IN_REGION);
	}
	
	public double getMaxMismatchRate() {
		return (Double) getOptions().valueOf(MAX_MISMATCH_RATE);
	}
	
	public boolean isPairedEnd() {
		return !getOptions().has(SINGLE_END);
	}
	
	public boolean isSkipAssembly() {
		return getOptions().has(SKIP_ASSEMBLY);
	}
	
	public boolean isSkipUnmappedAssemblyTrigger() {
		return getOptions().has(SKIP_UNMAPPED_ASSEMBLY_TRIGGER);
	}
	
	public boolean isKeepTmp() {
		return getOptions().has(KEEP_TMP);
	}
	
	public String getTmpDir() {
		return getOptions().has(TMP_DIR) ? (String) getOptions().valueOf(TMP_DIR) : null;
	}
	
	public String getInputVcf() {
		return getOptions().has(INPUT_VCF) ? (String) getOptions().valueOf(INPUT_VCF) : null;
	}
	
	public boolean useObservedIndels() {
		return !getOptions().has(SKIP_OBS_INDELS);
	}
	
	public boolean useSoftClippedReads() {
		return !getOptions().has(SKIP_SOFT_CLIP);
	}
	
	public boolean useConsensusSequence() {
		return getOptions().has(CONSENSUS_SEQ);
	}
	
	public boolean shouldUnsetDuplicates() {
		return getOptions().has(UNSET_DUPLICATES);
	}
	
	public int getMinBaseQuality() {
		return (Integer) getOptions().valueOf(MIN_BASE_QUALITY);
	}
	
	public double getMinReadCandidateFraction() {
		return (Double) getOptions().valueOf(MIN_READ_CANDIDATE_FRACTION);
	}
	
	public double getMinEdgeRatio() {
		return (Double) getOptions().valueOf(MIN_EDGE_RATIO);
	}
	
	public int getMaxAverageRegionDepth() {
		return (Integer) getOptions().valueOf(MAX_AVERAGE_REGION_DEPTH);
	}
	
	public int getMinimumMappingQuality() {
		return (Integer) getOptions().valueOf(MIN_MAPQ);
	}
	
	public boolean shouldSort() {
		return !getOptions().has(NO_SORT);
	}
	
	public boolean shouldCreateIndex() {
		return getOptions().has(INDEX);
	}
	
	public boolean shouldUseGkl() {
		return !getOptions().has(NO_GKL);
	}
	
	public int[] getSmithWatermanScoring() {
		String scoring = (String) getOptions().valueOf(SG_ALIGNER_SCORING);
		String[] fields = scoring.split(",");
		if (fields.length != 4) {
			Logger.error("4 values required for sga scoring");
			throw new IllegalArgumentException("4 values required for sga scoring");
		}
		
		int[] scores = new int[4];
		
		for (int i=0; i<4; i++) {
			scores[i] = Integer.parseInt(fields[i].trim());
		}
		
		return scores;
	}
	
	public int getWindowSize() {
		String win = (String) getOptions().valueOf(WINDOW_SIZE);
		String[] fields = win.split(",");
		if (fields.length != 2) {
			Logger.error("Please specify window size and overlap");
			throw new IllegalArgumentException("Window size and overlap must be specified");
		}
		
		return Integer.parseInt(fields[0].trim());
	}
	
	public int getWindowOverlap() {
		String win = (String) getOptions().valueOf(WINDOW_SIZE);
		String[] fields = win.split(",");
		if (fields.length != 2) {
			Logger.error("Please specify window size and overlap");
			throw new IllegalArgumentException("Window size and overlap must be specified");
		}
		
		return Integer.parseInt(fields[1].trim());
	}
	
	public int[] getContigAnchor() {
		String params = (String) getOptions().valueOf(CONTIG_ANCHOR);
		String[] fields = params.split(",");
		if (fields.length != 2) {
			Logger.error("2 values required for contig anchor");
		}
		
		int[] values = new int[2];
		
		for (int i=0; i<2; i++) {
			values[i] = Integer.parseInt(fields[i].trim());
		}	
		
		return values;
	}
	
	public int[] getSoftClipParams() {
		String params = (String) getOptions().valueOf(SOFT_CLIP);
		String[] fields = params.split(",");
		if (fields.length != 4) {
			Logger.error("4 values required for SW soft clip params");
		}
		
		int[] values = new int[4];
		
		for (int i=0; i<4; i++) {
			values[i] = Integer.parseInt(fields[i].trim());
		}	
		
		return values;
	}
	
	public String getLoggerLevel() {
		return (String) getOptions().valueOf(LOG_LEVEL);
	}
	
	public String getContigFile() {
		return (String) getOptions().valueOf(CONTIG_FILE);
	}
	
	public boolean isValid() {
		return isValid;
	}
	
	public int getMaxNodes() {
		return (Integer) getOptions().valueOf(MAX_NODES);
	}
	
	public int getMaxCachedReads() {
		return (Integer) getOptions().valueOf(MAX_CACHED_READS);
	}
	
	public int getCompressionLevel() {
		return (Integer) getOptions().valueOf(COMPRESSION_LEVEL);
	}
}
