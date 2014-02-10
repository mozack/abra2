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
	private static final String WORKING_DIR = "working";
	private static final String KMER_SIZE = "kmer";
	private static final String MIN_NODE_FREQUENCY = "mnf";
	private static final String MIN_UNALIGNED_NODE_FREQUENCY = "umnf";
	private static final String MIN_CONTIG_LENGTH = "mcl";
	private static final String MAX_POTENTIAL_CONTIGS = "mpc";
	private static final String MIN_CONTIG_MAPQ = "mc-mapq";
	private static final String NUM_THREADS = "threads";
	private static final String UNALIGNED_ASSEMBLY = "aur";
	private static final String MAX_UNALIGNED_READS = "mur";
	private static final String SINGLE_END = "single";
	private static final String RNA = "rna";
	private static final String RNA_OUTPUT = "rna-out";
	private static final String MIN_BASE_QUALITY = "mbq";
	private static final String MIN_READ_CANDIDATE_FRACTION = "rcf";
	
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required list of input sam or bam file(s) separated by comma").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM, "Required list of output sam or bam file(s) separated by comma").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE, "Genome reference location").withRequiredArg().ofType(String.class);
            parser.accepts(TARGET_REGIONS, "BED or GTF file containing target regions").withRequiredArg().ofType(String.class);
            parser.accepts(WORKING_DIR, "Working directory for intermediate output.  Must not already exist").withRequiredArg().ofType(String.class);
            parser.accepts(KMER_SIZE, "Assembly kmer size(delimit with commas if multiple sizes specified)").withRequiredArg().ofType(String.class);
            parser.accepts(MIN_NODE_FREQUENCY, "Assembly minimum node frequency").withRequiredArg().ofType(Integer.class).defaultsTo(2);
            parser.accepts(MIN_UNALIGNED_NODE_FREQUENCY, "Assembly minimum unaligned node frequency").withOptionalArg().ofType(Integer.class).defaultsTo(2);
            parser.accepts(MIN_CONTIG_LENGTH, "Assembly minimum contig length").withOptionalArg().ofType(Integer.class).defaultsTo(-1);
            parser.accepts(MAX_POTENTIAL_CONTIGS, "Maximum number of potential contigs for a region").withOptionalArg().ofType(Integer.class).defaultsTo(5000);
            parser.accepts(NUM_THREADS, "Number of threads").withRequiredArg().ofType(Integer.class).defaultsTo(2);
            parser.accepts(MIN_CONTIG_MAPQ, "Minimum contig mapping quality").withOptionalArg().ofType(Integer.class).defaultsTo(25);
            parser.accepts(UNALIGNED_ASSEMBLY, "Assemble unaligned reads (currently disabled).");
            parser.accepts(MAX_UNALIGNED_READS, "Maximum number of unaligned reads to assemble").withOptionalArg().ofType(Integer.class).defaultsTo(50000000);
            parser.accepts(SINGLE_END, "Input is single end");
            parser.accepts(RNA, "Input RNA sam or bam file (currently disabled)").withOptionalArg().ofType(String.class);
            parser.accepts(RNA_OUTPUT, "Output RNA sam or bam file (required if RNA input file specified)").withRequiredArg().ofType(String.class);
            parser.accepts(MIN_BASE_QUALITY, "Minimum base quality for inclusion in assembly").withOptionalArg().ofType(Integer.class).defaultsTo(20);
            parser.accepts(MIN_READ_CANDIDATE_FRACTION, "Minimum read candidate fraction for triggering assembly").withRequiredArg().ofType(Double.class).defaultsTo(.02);
    	}
    	
    	return parser;
	}

	@Override
	protected void validate() {
		isValid = true;
		
		if (!getOptions().hasArgument(INPUT_SAM)) {
			isValid = false;
			System.out.println("Missing required input SAM/BAM file");
		}

		if (!getOptions().hasArgument(OUTPUT_SAM)) {
			isValid = false;
			System.out.println("Missing required input SAM/BAM file");
		}
		
		if (getInputFiles().length != getOutputFiles().length) {
			System.out.println("Number of input files must equal number of output files");
		}
		
		if (!getOptions().hasArgument(REFERENCE)) {
			isValid = false;
			System.out.println("Missing required reference");
		}
		
		if (!getOptions().hasArgument(TARGET_REGIONS)) {
			isValid = false;
			System.out.println("Missing required target regions");
		}
		
		if (!getOptions().hasArgument(WORKING_DIR)) {
			isValid = false;
			System.out.println("Missing required working directory");
		}
		
		if ((getOptions().hasArgument(NUM_THREADS) && (Integer) getOptions().valueOf(NUM_THREADS) < 1)) {
			isValid = false;
			System.out.println("Num threads must be greater than zero.");
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
		return (String) getOptions().valueOf(TARGET_REGIONS);
	}
	
	public String getWorkingDir() {
		return (String) getOptions().valueOf(WORKING_DIR);
	}
	
	public int[] getKmerSizes() {
		String[] kmerStr = ((String) getOptions().valueOf(KMER_SIZE)).split(",");
		
		int[] kmers = new int[kmerStr.length];
		for (int i=0; i<kmerStr.length; i++) {
			kmers[i] = Integer.parseInt(kmerStr[i]);
		}
		
		return kmers;
	}
		
	public int getMinNodeFrequency() {
		return (Integer) getOptions().valueOf(MIN_NODE_FREQUENCY);
	}
	
	public int getMinUnalignedNodeFrequency() {
		return (Integer) getOptions().valueOf(MIN_UNALIGNED_NODE_FREQUENCY);
	}
	
	public int getMinContigLength() {
		return (Integer) getOptions().valueOf(MIN_CONTIG_LENGTH);
	}
	
	public int getMaxPotentialContigs() {
		return (Integer) getOptions().valueOf(MAX_POTENTIAL_CONTIGS);
	}
	
	public int getNumThreads() {
		return getOptions().hasArgument(NUM_THREADS) ? (Integer) getOptions().valueOf(NUM_THREADS) : 2;
	}
	
	public int getMinContigMapq() {
		return (Integer) getOptions().valueOf(MIN_CONTIG_MAPQ);
	}
	
	public boolean isSkipUnalignedAssembly() {
		return !getOptions().has(UNALIGNED_ASSEMBLY);
	}
	
	public int getMaxUnalignedReads() {
		return (Integer) getOptions().valueOf(MAX_UNALIGNED_READS);
	}
	
	public boolean isPairedEnd() {
		return !getOptions().has(SINGLE_END);
	}
	
	public String getRnaSam() {
		return (String) getOptions().valueOf(RNA);
	}
	
	public String getRnaSamOutput() {
		return (String) getOptions().valueOf(RNA_OUTPUT);
	}
	
	public int getMinBaseQuality() {
		return (Integer) getOptions().valueOf(MIN_BASE_QUALITY);
	}
	
	public double getMinReadCandidateFraction() {
		return (Double) getOptions().valueOf(MIN_READ_CANDIDATE_FRACTION);
	}
	
	public boolean isValid() {
		return isValid;
	}
}
