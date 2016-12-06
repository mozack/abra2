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
	private static final String MAX_POTENTIAL_CONTIGS = "mpc";
	private static final String MIN_MAPQ = "mapq";
	private static final String NUM_THREADS = "threads";
	private static final String SINGLE_END = "single";
	private static final String MIN_BASE_QUALITY = "mbq";
	private static final String MIN_READ_CANDIDATE_FRACTION = "rcf";
	private static final String MAX_AVERAGE_REGION_DEPTH = "mad";
	private static final String AVERAGE_DEPTH_CEILING = "adc";
	private static final String MIN_EDGE_RATIO = "mer";
	private static final String MAX_NODES = "maxn";
	private static final String SKIP_ASSEMBLY = "sa";
	private static final String SKIP_NON_ASSEMBLY = "sna";
	private static final String JUNCTIONS = "junctions";
	private static final String LOG_LEVEL = "log";
	private static final String CONTIG_FILE	 = "contigs";
	private static final String GTF_JUNCTIONS = "gtf";
	
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
            parser.accepts(MIN_NODE_FREQUENCY, "Assembly minimum node frequency").withRequiredArg().ofType(Integer.class).defaultsTo(2);
            parser.accepts(MIN_CONTIG_LENGTH, "Assembly minimum contig length").withOptionalArg().ofType(Integer.class).defaultsTo(-1);
            parser.accepts(MAX_POTENTIAL_CONTIGS, "Maximum number of potential contigs for a region").withOptionalArg().ofType(Integer.class).defaultsTo(5000);
            parser.accepts(NUM_THREADS, "Number of threads").withRequiredArg().ofType(Integer.class).defaultsTo(4);
            parser.accepts(MIN_MAPQ, "Minimum mapping quality for a read to be used in assembly and be eligible for realignment").withOptionalArg().ofType(Integer.class).defaultsTo(20);
            parser.accepts(SINGLE_END, "Input is single end");
            parser.accepts(MIN_BASE_QUALITY, "Minimum base quality for inclusion in assembly.  This value is compared against the sum of base qualities per kmer position").withOptionalArg().ofType(Integer.class).defaultsTo(60);
            parser.accepts(MIN_READ_CANDIDATE_FRACTION, "Minimum read candidate fraction for triggering assembly").withRequiredArg().ofType(Double.class).defaultsTo(.01);
            parser.accepts(MAX_AVERAGE_REGION_DEPTH, "Regions with average depth exceeding this value will be downsampled").withRequiredArg().ofType(Integer.class).defaultsTo(250);
            parser.accepts(AVERAGE_DEPTH_CEILING, "Skip regions with average depth greater than this value").withOptionalArg().ofType(Integer.class).defaultsTo(100000);
            parser.accepts(MIN_EDGE_RATIO, "Min edge pruning ratio.  Default value is appropriate for relatively sensitive somatic cases.  May be increased for improved speed in germline only cases.").withRequiredArg().ofType(Double.class).defaultsTo(.02);
            parser.accepts(MAX_NODES, "Maximum pre-pruned nodes in regional assembly").withOptionalArg().ofType(Integer.class).defaultsTo(9000);
            parser.accepts(SKIP_ASSEMBLY, "Skip assembly");
            parser.accepts(SKIP_NON_ASSEMBLY, "Skip non-assembly contig generation");
            parser.accepts(JUNCTIONS, "Splice junctions definition file").withRequiredArg().ofType(String.class);
            parser.accepts(LOG_LEVEL, "Logging level (trace,debug,info,warn,error)").withRequiredArg().ofType(String.class).defaultsTo("info");
            parser.accepts(CONTIG_FILE, "Optional file to which assembled contigs are written").withRequiredArg().ofType(String.class);
            parser.accepts(GTF_JUNCTIONS, "GTF file defining exons and transcripts").withRequiredArg().ofType(String.class);
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
	
	public int getMaxPotentialContigs() {
		return (Integer) getOptions().valueOf(MAX_POTENTIAL_CONTIGS);
	}
	
	public int getNumThreads() {
		return getOptions().hasArgument(NUM_THREADS) ? (Integer) getOptions().valueOf(NUM_THREADS) : 4;
	}
	
	public boolean isPairedEnd() {
		return !getOptions().has(SINGLE_END);
	}
	
	public boolean isSkipAssembly() {
		return getOptions().has(SKIP_ASSEMBLY);
	}
	
	public boolean isSkipNonAssembly() {
		return getOptions().has(SKIP_NON_ASSEMBLY);
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
	
	public int getAverageDepthCeiling() {
		return (Integer) getOptions().valueOf(AVERAGE_DEPTH_CEILING);
	}
	
	public int getMinimumMappingQuality() {
		return (Integer) getOptions().valueOf(MIN_MAPQ);
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
}
