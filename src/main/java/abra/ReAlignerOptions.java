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
	private static final String INPUT_SAM2 = "in2";
	private static final String INPUT_SAM3 = "in3";
	private static final String OUTPUT_SAM = "out";
	private static final String OUTPUT_SAM2 = "out2";
	private static final String OUTPUT_SAM3 = "out3";
	private static final String REFERENCE = "ref";
	private static final String TARGET_REGIONS = "targets";
	private static final String WORKING_DIR = "working";
	private static final String KMER_SIZE = "kmer";
	private static final String MIN_NODE_FREQUENCY = "mnf";
	private static final String MIN_UNALIGNED_NODE_FREQUENCY = "umnf";
	private static final String MIN_CONTIG_LENGTH = "mcl";
	private static final String MAX_POTENTIAL_CONTIGS = "mpc";
	private static final String MIN_CONTIG_RATIO = "mcr";
	private static final String MIN_CONTIG_MAPQ = "mc-mapq";
	private static final String NUM_THREADS = "threads";
	private static final String SKIP_UNALIGNED_ASSEMBLY = "no-unalign";
	private static final String MAX_UNALIGNED_READS = "mur";
	private static final String PAIRED_END = "paired";
	private static final String RNA = "rna";
	private static final String RNA_OUTPUT = "rna-out";
	private static final String MIN_BASE_QUALITY = "mbq";
	private static final String FAVOR_GAP_EXTENSIONS = "fgap";
	private static final String FILTER_SNP_CLUSTERS = "fsc";
	private static final String PAD_REGIONS = "pad";
	
	private OptionParser parser;
	private boolean isValid;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(INPUT_SAM, "Required input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(INPUT_SAM2, "Optional input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(INPUT_SAM3, "Optional input sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM, "Required output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM2, "Optional output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(OUTPUT_SAM3, "Optional output sam or bam file").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE, "Genome reference location").withRequiredArg().ofType(String.class);
            parser.accepts(TARGET_REGIONS, "GTF containing target regions").withRequiredArg().ofType(String.class);
            parser.accepts(WORKING_DIR, "Working directory for intermediate output").withRequiredArg().ofType(String.class);
            parser.accepts(KMER_SIZE, "Assembly kmer size(delimit by commas if more than 1").withRequiredArg().ofType(String.class);
            parser.accepts(MIN_NODE_FREQUENCY, "Assembly minimum node frequency").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_UNALIGNED_NODE_FREQUENCY, "Assembly minimum unaligned node frequency").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_CONTIG_LENGTH, "Assembly minimum contig length").withRequiredArg().ofType(Integer.class);
            parser.accepts(MAX_POTENTIAL_CONTIGS, "Maximum number of potential contigs for a region").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_CONTIG_RATIO, "Minimum contig length as percentage of observed region length").withRequiredArg().ofType(Double.class);
            parser.accepts(NUM_THREADS, "Number of threads (default: 2)").withRequiredArg().ofType(Integer.class);
            parser.accepts(MIN_CONTIG_MAPQ, "Minimum contig mapping quality").withRequiredArg().ofType(Integer.class);
            parser.accepts(SKIP_UNALIGNED_ASSEMBLY, "Skip assembly of reads that did not initially align.");
            parser.accepts(MAX_UNALIGNED_READS, "Maximum number of unaligned reads to assemble").withRequiredArg().ofType(Integer.class);
            parser.accepts(PAIRED_END, "Paired end");
            parser.accepts(RNA, "Input RNA sam or bam file (optional)").withRequiredArg().ofType(String.class);
            parser.accepts(RNA_OUTPUT, "Output RNA sam or bam file (required if RNA input file specified)").withRequiredArg().ofType(String.class);
            parser.accepts(MIN_BASE_QUALITY, "Minimum base quality for inclusion in assembly").withRequiredArg().ofType(Integer.class);
            parser.accepts(FAVOR_GAP_EXTENSIONS, "Reduce penalty for gap extensions");
            parser.accepts(FILTER_SNP_CLUSTERS, "Filter high numbers of nearby mismatches");
            parser.accepts(PAD_REGIONS, "Expand target regions by read length");
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
	
	public String getInputFile() {
		return (String) getOptions().valueOf(INPUT_SAM);
	}
	
	public String getOutputFile() {
		return (String) getOptions().valueOf(OUTPUT_SAM);
	}
	
	public String getInputFile2() {
		return (String) getOptions().valueOf(INPUT_SAM2);
	}
	
	public String getOutputFile2() {
		return (String) getOptions().valueOf(OUTPUT_SAM2);
	}
	
	public String getInputFile3() {
		return (String) getOptions().valueOf(INPUT_SAM3);
	}
	
	public String getOutputFile3() {
		return (String) getOptions().valueOf(OUTPUT_SAM3);
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
	
	public double getMinContigRatio() {
		return (Double) getOptions().valueOf(MIN_CONTIG_RATIO);
	}
	
	public int getNumThreads() {
		return getOptions().hasArgument(NUM_THREADS) ? (Integer) getOptions().valueOf(NUM_THREADS) : 2;
	}
	
	public int getMinContigMapq() {
		return (Integer) getOptions().valueOf(MIN_CONTIG_MAPQ);
	}
	
	public boolean isSkipUnalignedAssembly() {
		return getOptions().has(SKIP_UNALIGNED_ASSEMBLY);
	}
	
	public int getMaxUnalignedReads() {
		return (Integer) getOptions().valueOf(MAX_UNALIGNED_READS);
	}
	
	public boolean isPairedEnd() {
		return getOptions().has(PAIRED_END);
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
	
	public boolean isGapExtensionFavored() {
		return getOptions().has(FAVOR_GAP_EXTENSIONS);
	}
	
	public boolean isFilterSnpClusters() {
		return getOptions().has(FILTER_SNP_CLUSTERS);
	}
	
	public boolean isPadRegions() {
		return getOptions().has(PAD_REGIONS);
	}
	
	public boolean isValid() {
		return isValid;
	}
}
