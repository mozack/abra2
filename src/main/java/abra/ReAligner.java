/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import static abra.Logger.log;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import picard.sam.BuildBamIndex;
import picard.sam.SortSam;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

/**
 * ABRA's main entry point
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAligner {

	private static final int DEFAULT_MAX_UNALIGNED_READS = 1000000;
	
	public static final int MAX_REGION_LENGTH = 400;
	private static final int MIN_REGION_REMAINDER = 200;
	public static final int REGION_OVERLAP = 200;

	private static final int MAX_POTENTIAL_UNALIGNED_CONTIGS = 2000000;
	
	// Minimum sequence length recommended for use with bwa mem
	private static final int MIN_CONTIG_LENGTH = 70;
	
	private SAMFileHeader[] samHeaders;
	
	private List<Feature> regions;

	private String regionsBed;

	private String tempDir;
	
	private String unalignedRegionSam;

	private String reference;
	
	private String bwaIndex;
	
	private int minContigMapq;

	private AssemblerSettings assemblerSettings;
	
	private int numThreads;
	
	private int maxUnalignedReads = DEFAULT_MAX_UNALIGNED_READS;
	
	private boolean shouldReprocessUnaligned = true;
	
	private boolean isOutputIntermediateBam = false;
	
	private String structuralVariantFile;
	private String localRepeatFile;
	private BufferedWriter localRepeatWriter;
	
	private String[] inputSams;
	private SAMFileWriter[] writers;
	private String[] tempDirs;
	
	private int readLength = -1;
	private int maxMapq = -1;
	private int minInsertLength = Integer.MAX_VALUE;
	private int maxInsertLength = -1;
	
	private boolean isPairedEnd = false;
	
	private String rnaSam = null;
	private String rnaOutputSam = null;
	private SAMFileHeader rnaHeader = null;
	private int rnaReadLength = -1;
	
	private BufferedWriter contigWriter;
	private BufferedWriter svContigWriter;
	
	private CompareToReference2 c2r;
	
	private ReadAdjuster readAdjuster;
	
	private ThreadManager threadManager;
	
	private boolean hasContigs = false;
	
	private int minMappingQuality;
	
	private boolean isDebug;
	
	// If true, the input target file specifies kmer values
	private boolean hasPresetKmers = false;
	
	public static final int COMPRESSION_LEVEL = 1;
	
	public void reAlign(String[] inputFiles, String[] outputFiles) throws Exception {
		
		this.inputSams = inputFiles;
		
		logStartupInfo(outputFiles);
				
		init();
		
		c2r = new CompareToReference2();
		c2r.init(this.reference);

		log("Reading Input SAM Header and identifying read length");
		getSamHeaderAndReadLength();
		
		log("Read length: " + readLength);
		
		log("Loading target regions");
		loadRegions();
		
		readAdjuster = new ReadAdjuster(isPairedEnd, this.maxMapq, c2r, minInsertLength, maxInsertLength);
		
		if (shouldReprocessUnaligned) {
//			processUnaligned();
		}
		
		Clock clock = new Clock("Assembly");
		clock.start();
		
		String contigFasta = tempDir + "/" + "all_contigs.fasta";
		contigWriter = new BufferedWriter(new FileWriter(contigFasta, false));
		
		String svContigFasta = tempDir + "/" + "sv_contigs.fasta";
		svContigWriter = new BufferedWriter(new FileWriter(svContigFasta, false));
		
		tempDirs = new String[inputSams.length];
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
//		writerFactory.setUseAsyncIo(true);
//		writerFactory.setAsyncOutputBufferSize(500000);
		
		writerFactory.setUseAsyncIo(false);
		
		writers = new SAMFileWriter[inputSams.length];
		
		for (int i=0; i<inputSams.length; i++) {
			// init temp dir
			String temp = tempDir + "/temp" + (i+1);
			mkdir(temp);
			tempDirs[i] = temp;

			// init BAM writer
			writers[i] = writerFactory.makeBAMWriter(
					samHeaders[i], false, new File(outputFiles[i]), COMPRESSION_LEVEL);
		}
		
		// Start pre-processing reads on separate thread for each input file.
		// This happens in parallel with assembly, provided there are enough threads.
		for (int i=0; i<inputSams.length; i++) {
			preProcessReads(inputSams[i], tempDirs[i], writers[i]);
		}
		
		log("Iterating over regions");
		
		for (Feature region : regions) {
			spawnRegionThread(region, null);
		}
		
		
		log("Waiting for all threads to complete");
		threadManager.waitForAllThreadsToComplete();
		
		contigWriter.close();
		svContigWriter.close();
		
		if (localRepeatWriter != null) {
			localRepeatWriter.close();
		}
		
		clock.stopAndPrint();
		
		String cleanContigsFasta = null;
		
		if (hasContigs) {
			clock = new Clock("Align and clean contigs");
			clock.start();
			cleanContigsFasta = alignAndCleanContigs(contigFasta, tempDir, true);
			clock.stopAndPrint();
		}
				
		if (cleanContigsFasta != null) {		
			clock = new Clock("Align to contigs");
			clock.start();
			
			String[] alignedSams = alignReads(cleanContigsFasta, c2r);
			
			clock.stopAndPrint();
			
			for (SAMFileWriter writer : this.writers) {
				writer.close();
			}
			
			if (rnaSam != null) {
				processRna();
			}
			
		} else {
			log("WARNING!  No contigs assembled.  Just making a copy of input converting to/from SAM/BAM as appropriate.");
			for (int i=0; i<inputFiles.length; i++) {
				copySam(inputFiles[i], outputFiles[i]);	
			}
		}
		
		if (this.assemblerSettings.searchForStructuralVariation() && this.isPairedEnd) {
			clock = new Clock("Structural Variant search");
			clock.start();
			new SVEvaluator().evaluateAndOutput(svContigFasta, this, tempDir, readLength, inputFiles, tempDirs, samHeaders, structuralVariantFile);
			clock.stopAndPrint();
		}
		
		System.out.println("Done.");
	}
	
	private void preProcessReads(String inputSam, String tempDir, SAMFileWriter writer) throws InterruptedException {
		PreprocessReadsRunnable thread = new PreprocessReadsRunnable(threadManager, this,
				inputSam, this.getTempReadFile(tempDir), c2r, writer);

		threadManager.spawnThread(thread);
	}
	
	private void processRna() {
		/*
		if (rnaSam != null) {
			String rnaTemp = tempDir + "/rna";
			mkdir(rnaTemp);
			getRnaSamHeaderAndReadLength(rnaSam);
			rnaHeader.setSortOrder(SortOrder.coordinate);
			
			clock = new Clock("RNA - Sam2Fastq and Align");
			clock.start();
			log("Aligning RNA to contigs");
			String alignedToContigRna = alignReads(rnaTemp, rnaSam, cleanContigsFasta, c2r, rnaOutputSam);
			clock.stopAndPrint();
			
			clock = new Clock("Adjust reads");
			clock.start();
			log("Adjusting RNA reads");
			adjustReads(alignedToContigRna, rnaOutputSam, true, c2r);
			clock.stopAndPrint();
		}
		*/
	}
	
	private void logStartupInfo(String[] outputFiles) {
		
		int ctr = 0;
		for (String input : inputSams) {
			System.out.println("input" + ctr + ": " + input);
		}

		ctr = 0;
		for (String output : outputFiles) {
			System.out.println("output" + ctr + ": " + output);
		}
		
		System.out.println("regions: " + regionsBed);
		System.out.println("reference: " + reference);
		System.out.println("bwa index: " + bwaIndex);
		System.out.println("working dir: " + tempDir);
		System.out.println("num threads: " + numThreads);
		System.out.println("max unaligned reads: " + maxUnalignedReads);
		System.out.println(assemblerSettings.getDescription());
		System.out.println("rna: " + rnaSam);
		System.out.println("rna output: " + rnaOutputSam);
		System.out.println("paired end: " + isPairedEnd);
		System.out.println("use intermediate bam: " + isOutputIntermediateBam);
		
		String javaVersion = System.getProperty("java.version");
		System.out.println("Java version: " + javaVersion);
		if (javaVersion.startsWith("1.6") || javaVersion.startsWith("1.5") || javaVersion.startsWith("1.4")) {
			throw new RuntimeException("Please upgrade to Java 7 or later to run ABRA.");
		}
		
		try {
			InetAddress localhost = java.net.InetAddress.getLocalHost();
			String hostname = localhost.getHostName();
			System.out.println("hostname: " + hostname);
		} catch (Throwable t) {
			System.out.println("Error getting hostname: " + t.getMessage());
		}
	}
	
	/*
	private void processUnaligned() throws IOException, InterruptedException {
		Clock clock = new Clock("Process unaligned");
		clock.start();
		log("Assembling unaligned reads");
		
		String unalignedSam = tempDir + "/unaligned.bam";
		unalignedSam = getUnalignedReads(unalignedSam);
		
//		String unalignedSam = tempDir + "/" + "unaligned_to_contig.bam";
		
		String unalignedDir = tempDir + "/unaligned";
		String unalignedContigFasta = unalignedDir + "/unaligned_contigs.fasta";
		unalignedRegionSam = unalignedDir + "/unaligned_region.bam";
		String sortedUnalignedRegion = unalignedDir + "/sorted_unaligned_region.bam";
		
		Assembler assem = newUnalignedAssembler(1);
		List<String> unalignedSamList = new ArrayList<String>();
		unalignedSamList.add(unalignedSam);
		List<String> unalignedAssemblies = assem.assembleContigs(unalignedSamList, unalignedContigFasta, tempDir, null, "unaligned", false, this, null);
		
		boolean hasContigs = unalignedAssemblies.size() > 0;

		// Make eligible for GC
		assem = null;
					
		if (hasContigs) {
			unalignedContigFasta = unalignedAssemblies.get(0);
			String unalignedCleanContigsFasta = alignAndCleanContigs(unalignedContigFasta, unalignedDir, false);
			if (unalignedCleanContigsFasta != null) {
				// Build contig fasta index
				log("Indexing contigs from unaligned region");
				Aligner contigAligner = new Aligner(unalignedCleanContigsFasta, numThreads);
				contigAligner.index();
				log("Done Indexing contigs from unaligned region");
				String alignedToContigSam = unalignedDir + "/" + "align_to_contig.bam";
				alignReads(unalignedDir, unalignedSam, unalignedCleanContigsFasta, null, null, alignedToContigSam);
				String alignedToContigBam = alignedToContigSam;
				
				log("Adjusting unaligned reads");
				SAMFileWriter unalignedWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
						samHeader, false, new File(unalignedRegionSam));
				ReadAdjuster unalignedReadAdjuster = new ReadAdjuster(isPairedEnd, this.maxMapq, null, minInsertLength, maxInsertLength);					
				unalignedReadAdjuster.adjustReads(alignedToContigBam, unalignedWriter, false, unalignedDir, samHeader);
				unalignedWriter.close();
				
				sortBam(unalignedRegionSam, sortedUnalignedRegion, "coordinate");
				unalignedRegionSam = sortedUnalignedRegion;
				
				indexBam(unalignedRegionSam);
			} else {
				shouldReprocessUnaligned = false;
			}
		} else {
			shouldReprocessUnaligned = false;
		}
		clock.stopAndPrint();

	}
	*/
	
	private void copySam(String input, String output) {
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(output));
		
		for (SAMRecord read : reader) {
			writer.addAlignment(read);
		}
		
		reader.close();
		writer.close();
	}
	
	private String[] alignReads(String cleanContigsFasta, CompareToReference2 c2r) throws IOException, InterruptedException {
		
		// Build contig fasta index
		log("Indexing contigs");
		Aligner contigAligner = new Aligner(cleanContigsFasta, numThreads);
		contigAligner.index();
		log("Contig indexing done");
		
		String[] alignedToContigsSams = new String[inputSams.length];
		
		for (int i=0; i<inputSams.length; i++) {
			alignedToContigsSams[i] = tempDirs[i] + "/" + "align_to_contig.sam";
			alignReads(tempDirs[i], inputSams[i], cleanContigsFasta, c2r, writers[i], alignedToContigsSams[i], samHeaders[i]);			
		}
				
		return alignedToContigsSams;
	}
	
	public static int getNumIndelBases2(SAMRecord read) {
		int numIndelBases = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if (element.getOperator() == CigarOperator.D) {
				numIndelBases += 1;
			} else if (element.getOperator() == CigarOperator.I) {
				numIndelBases += element.getLength();
			}
		}
		
		return numIndelBases;
	}

	private void discardMisalignedContigs(String inputSam, String outputSam) {
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter outputReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeaders[0], true, new File(outputSam));

		for (SAMRecord contig : reader) {
			String[] fields = contig.getReadName().split("_");
			
			String regionChromosome = "";
			
			// Loop through fields in case the chromosome name contains
			// an underscore.
			for (int i=0; i<fields.length-3; i++) {
				regionChromosome += fields[i];
				if (i+1 < fields.length-3) {
					regionChromosome += "_";
				}
			}
			
			int regionStart = Integer.parseInt(fields[fields.length-3]) - 1000;
			int regionStop = Integer.parseInt(fields[fields.length-2]) + 1000;
						
			if ((contig.getReferenceName().equals(regionChromosome)) &&
				(contig.getAlignmentStart() >= regionStart) &&
				(contig.getAlignmentEnd() <= regionStop)) {
			
				outputReadsBam.addAlignment(contig);
			}
		}
		
		outputReadsBam.close();
		reader.close();
	}
	
	void alignStructuralVariantCandidates(String svContigFasta, String svContigsSam) throws InterruptedException, IOException {
		Aligner aligner = new Aligner(reference, numThreads);
		aligner.align(svContigFasta, svContigsSam, false);
	}
	
	private String alignAndCleanContigs(String contigFasta, String tempDir, boolean isTightAlignment) throws InterruptedException, IOException {
		log("Aligning contigs");
		Aligner aligner = new Aligner(bwaIndex, numThreads);
		String contigsSam = tempDir + "/" + "all_contigs.sam";
		aligner.align(contigFasta, contigsSam, false);
		
		if (isTightAlignment) {
			log("Discarding contigs aligned outside of region");
			String allInRegionSam = tempDir + "/" + "all_contigs_in_region.sam";
			discardMisalignedContigs(contigsSam, allInRegionSam);
			contigsSam = allInRegionSam;
		}
		
		log("Processing chimeric reads");
		CombineChimera3 cc = new CombineChimera3();
		String contigsWithChim = tempDir + "/" + "all_contigs_chim.bam";
		int slack = this.readLength / 3;
		cc.combine(contigsSam, contigsWithChim, isTightAlignment ? slack : 0, c2r);
		
		if (isTightAlignment) {
			// Chop and clop...
			log("Sorting and indexing for chopper clopper.");
			String contigsWithChimSorted = tempDir + "/" + "all_contigs_chim_sorted.bam";
			sortBam(contigsWithChim, contigsWithChimSorted, "coordinate");
			indexBam(contigsWithChimSorted);
			
			log("Chopper clopper start.");
			ContigChopper chopper = new ContigChopper();
			chopper.setC2R(c2r);
			chopper.setReadLength(this.readLength);
			
			String contigsWithChimChopped = tempDir + "/" + "all_contigs_chim_chopped.bam";
			chopper.chopClopDrop(this.regions, contigsWithChimSorted, contigsWithChimChopped);
			
			log("Chopper clopper done.");
			contigsWithChim = contigsWithChimChopped;
			
			chopper = null;
		}
		
		log("Cleaning contigs");
		String cleanContigsFasta = tempDir + "/" + "clean_contigs.fasta";
		boolean hasCleanContigs = cleanAndOutputContigs(contigsWithChim, cleanContigsFasta, isTightAlignment);
		
		return hasCleanContigs ? cleanContigsFasta : null;
	}
	
	String alignReads(String tempDir, String inputSam, String cleanContigsFasta,
			CompareToReference2 c2r, SAMFileWriter finalOutputSam, String alignedToContigSam,
			SAMFileHeader header) throws InterruptedException, IOException {
		log("Aligning original reads to contigs");
		alignToContigs(tempDir, alignedToContigSam, cleanContigsFasta, finalOutputSam, header);
		return alignedToContigSam;
	}
	
	private void indexBam(String bam) {
		String[] args = new String[] { 
				"INPUT=" + bam,
				"VALIDATION_STRINGENCY=SILENT"
				};
		
		int ret = new BuildBamIndex().instanceMain(args);
		if (ret != 0) {
			throw new RuntimeException("BuildBamIndex failed");
		}
	}	
		
	void sortBam(String input, String output, String sortOrder) {
		String[] args = new String[] { 
				"INPUT=" + input, 
				"OUTPUT=" + output, 
				"VALIDATION_STRINGENCY=SILENT",
				"SORT_ORDER=" + sortOrder,
				"TMP_DIR=" + this.tempDir + "/sorttmp"
				};
		
		int ret = new SortSam().instanceMain(args);
		if (ret != 0) {
			throw new RuntimeException("SortSam failed");
		}
	}
	
	private void spawnRegionThread(Feature region, String inputSam) throws InterruptedException {
		ReAlignerRunnable thread = new ReAlignerRunnable(threadManager, this, region);
		threadManager.spawnThread(thread);
	}
		
	private boolean shouldIncludeInUnalignedPile(SAMRecord read) {
		boolean shouldInclude = false;
		
		if (!read.getReadFailsVendorQualityCheckFlag()) {
			if (read.getReadUnmappedFlag()) {
				shouldInclude = true;
			}
			// For Stampy, if Cigar length > 4 and read is not ambiguous (mapq >= 4)
			else if ((read.getCigarLength() > 4) && read.getMappingQuality() >= 4) {
				shouldInclude = true;
			}
		}
		
		return shouldInclude;
	}
	
	private synchronized void appendContigs(String contigs) throws IOException {
		contigWriter.write(contigs);
		hasContigs = true;
	}
	
	public void processRegion(Feature region) throws Exception {
		if (isDebug) {
			log("Processing region: " + region.getDescriptor());
		}
		
		try {
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			
			List<String> bams = new ArrayList<String>(Arrays.asList(this.inputSams));
			if (shouldReprocessUnaligned) {
				bams.add(unalignedRegionSam);
			}
			
			// Assemble contigs
			if (region.getKmer() > this.readLength-15) {
				System.out.println("Skipping assembly of region: " + region.getDescriptor() + " - " + region.getKmer());
			} else {
				NativeAssembler assem = (NativeAssembler) newAssembler(region);
				List<Feature> regions = new ArrayList<Feature>();
				regions.add(region);
				String contigs = assem.assembleContigs(bams, contigsFasta, tempDir, regions, region.getDescriptor(), true, this, c2r);
				if (!contigs.equals("<ERROR>") && !contigs.equals("<REPEAT>") && !contigs.isEmpty()) {
					
					appendContigs(contigs);
				
					List<BreakpointCandidate> svCandidates = assem.getSvCandidateRegions();
					for (BreakpointCandidate svCandidate : svCandidates) {
						
						if (isDebug) {
							System.out.println("SV: " + region.getDescriptor() + "-->" + svCandidate.getRegion().getDescriptor());
						}
						List<Feature> svRegions = new ArrayList<Feature>();
						svRegions.add(region);
						Feature svCandidateRegion = new Feature(svCandidate.getRegion().getSeqname(), svCandidate.getRegion().getStart(), 
								Math.min(svCandidate.getRegion().getEnd(), c2r.getReferenceLength(svCandidate.getRegion().getSeqname())-1));
						
						svRegions.add(svCandidateRegion);
						
						NativeAssembler svAssem = (NativeAssembler) newAssembler(region);
						String svContigs = svAssem.assembleContigs(bams, contigsFasta, tempDir, svRegions, region.getDescriptor() + "__" + svCandidate.getRegion().getDescriptor() + "_" + svCandidate.getSpanningReadPairCount(), true, this, c2r);
						
						if (!svContigs.equals("<ERROR>") && !svContigs.equals("<REPEAT>") && !svContigs.isEmpty()) {
							svContigWriter.write(svContigs);
						}
					}
				}
				
				if (assem.isCycleExceedingThresholdDetected() && (bams.size() > 1) && this.localRepeatWriter != null) {
					System.out.println("Attempting cycle detection for: " + region.getDescriptor());
					
					// Assemble each region separately looking for cycles
					List<String> cycleStatus = new ArrayList<String>();
					
					int kmer = readLength / 2;
					if (kmer % 2 == 1) {
						kmer -= 1;
					}
					kmer = Math.min(kmer, NativeAssembler.CYCLE_KMER_LENGTH_THRESHOLD);
					kmer = Math.max(kmer, region.getKmer());
					
					for (String bam : bams) {
						List<String> bamInput = new ArrayList<String>();
						bamInput.add(bam);
						NativeAssembler cycleAssem = (NativeAssembler) newAssembler(region);
												
						cycleAssem.setKmer(new int[] { kmer });
						cycleAssem.setShouldSearchForSv(false);
						
						// Double pruning thresholds for repeat discovery
						cycleAssem.setMinBaseQuality(assemblerSettings.getMinBaseQuality() * 2);
						cycleAssem.setMinKmerFrequency(assemblerSettings.getMinNodeFrequncy() * 2);
						
						cycleAssem.setMinEdgeRatio(assemblerSettings.getMinEdgeRatio());
						
						String cycleContigs = cycleAssem.assembleContigs(bamInput, contigsFasta, tempDir, regions, region.getDescriptor(), true, this, c2r);
						
						if (!cycleContigs.equals("<ERROR>") && !cycleContigs.equals("<REPEAT>")) {
							cycleContigs = ".";
						}
						
						cycleStatus.add(cycleContigs);
					}
					
					StringBuffer buf = new StringBuffer("Cycle detection result");
					for (String status : cycleStatus) {
						buf.append(status + "\t");
					}
					
					System.out.println("Cycle detection for region: " + region + ".  Result: [" + buf.toString() + "]");
					
					if (isAnyElementDifferent(cycleStatus)) {
						StringBuffer out = new StringBuffer();
						out.append(region.getDescriptor() + "\t");
						
						for (String status : cycleStatus) {
							out.append(status + "\t");							
						}
						
						out.append(kmer);
						out.append('\n');
						
						synchronized (localRepeatWriter) {
							localRepeatWriter.write(out.toString());
						}
					}
				}
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}
		
	private boolean isAnyElementDifferent(List<String> elems) {
		String last = null;
		
		for (String elem : elems) {
			if (last != null && !elem.equals(last)) {
				return true;
			}
			
			last = elem;
		}
		
		return false;
	}
	
	static List<Feature> getRegions(String regionsBed, int readLength, boolean hasPresetKmers) throws IOException {
		RegionLoader loader = new RegionLoader();
		List<Feature> regions = loader.load(regionsBed, hasPresetKmers);
		if (regions.size() > 0 && (regions.get(0).getKmer() == 0)) {
			regions = RegionLoader.collapseRegions(regions, readLength);
			regions = splitRegions(regions);
		}

		return regions;
	}
		
	private void loadRegions() throws IOException {
		this.regions = getRegions(regionsBed, readLength, hasPresetKmers);
		
		System.out.println("Num regions: " + regions.size());
		if (isDebug) {
			for (Feature region : regions) {
				System.out.println(region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd() + "\t" + region.getKmer());
			}
		}
	}

	public void setRegionsBed(String bedFile) {
		this.regionsBed = bedFile;
	}

	private void getSamHeaderAndReadLength() {
		
		log("Identifying header and determining read length");
		
		this.samHeaders = new SAMFileHeader[this.inputSams.length];
		
		for (int i=0; i<this.inputSams.length; i++) {
		
			SAMFileReader reader = new SAMFileReader(new File(inputSams[i]));
			try {
				reader.setValidationStringency(ValidationStringency.SILENT);
		
				samHeaders[i] = reader.getFileHeader();
				samHeaders[i].setSortOrder(SAMFileHeader.SortOrder.unsorted);
				
				Iterator<SAMRecord> iter = reader.iterator();
				
				int cnt = 0;
				while ((iter.hasNext()) && (cnt < 1000000)) {
					SAMRecord read = iter.next();
					this.readLength = Math.max(this.readLength, read.getReadLength());
					this.maxMapq = Math.max(this.maxMapq, read.getMappingQuality());
					
					// Assumes aligner sets proper pair flag correctly
					if ((isPairedEnd) && (read.getReadPairedFlag()) && (read.getProperPairFlag())) {
						this.minInsertLength = Math.min(this.minInsertLength, Math.abs(read.getInferredInsertSize()));
						this.maxInsertLength = Math.max(this.maxInsertLength, Math.abs(read.getInferredInsertSize()));
					}
				}
				
				// Allow some fudge in insert length
				minInsertLength = Math.max(minInsertLength - 2*readLength, 0);
				maxInsertLength = maxInsertLength + 2*readLength;
				
			} finally {
				reader.close();
			}
		}
		
		System.out.println("Min insert length: " + minInsertLength);
		System.out.println("Max insert length: " + maxInsertLength);
				
		log("Max read length is: " + readLength);
		if (assemblerSettings.getMinContigLength() < 1) {
			assemblerSettings.setMinContigLength(Math.max(readLength+1, MIN_CONTIG_LENGTH));
		}
		log("Min contig length: " + assemblerSettings.getMinContigLength());
	}
	
	//TODO: Dedup with getSamHeaderAndReadLength
	private void getRnaSamHeaderAndReadLength(String inputSam) {
		
		log("Identifying RNA header and determining read length");
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		try {
			reader.setValidationStringency(ValidationStringency.SILENT);
	
			rnaHeader = reader.getFileHeader();
			rnaHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
			
			Iterator<SAMRecord> iter = reader.iterator();
			
			int cnt = 0;
			while ((iter.hasNext()) && (cnt < 1000000)) {
				SAMRecord read = iter.next();
				this.rnaReadLength = Math.max(this.rnaReadLength, read.getReadLength());
			}
		} finally {
			reader.close();
		}
		
		log("Max RNA read length is: " + rnaReadLength);
	}
	
	void sam2Fastq(String bam, String intermediateOutput, CompareToReference2 c2r, SAMFileWriter finalOutputSam) throws IOException {
		log("Preprocessing: " + bam);
		Sam2Fastq sam2Fastq = new Sam2Fastq();
		sam2Fastq.convert(bam, intermediateOutput, c2r, finalOutputSam, isPairedEnd, regions, minMappingQuality, isOutputIntermediateBam);
		log("Done Preprocessing: " + bam);
	}
			
	private boolean cleanAndOutputContigs(String contigsSam, String cleanContigsFasta, boolean shouldRemoveSoftClips) throws IOException {
		
		boolean hasCleanContigs = false;
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(cleanContigsFasta, false));
		
		SAMFileReader contigReader = new SAMFileReader(new File(contigsSam));
		contigReader.setValidationStringency(ValidationStringency.SILENT);
		
		int contigCount = 0;
		
		for (SAMRecord contigRead : contigReader) {
			if (contigRead.getMappingQuality() >= this.minContigMapq) {
				
				SAMRecordUtils.removeSoftClips(contigRead);
				
				String bases = contigRead.getReadString();
				
				//TODO: Why would cigar length be zero here?
				if ((bases.length() >= assemblerSettings.getMinContigLength()) &&
					(contigRead.getCigarLength() > 0)) {
					
					CigarElement first = contigRead.getCigar().getCigarElement(0);
					CigarElement last = contigRead.getCigar().getCigarElement(contigRead.getCigarLength()-1);
					
					if ((first.getOperator() == CigarOperator.M) &&
						(last.getOperator() == CigarOperator.M)) {
						
						String prefix = "";
						String suffix = "";

						if (!contigRead.getReferenceName().startsWith("uc0")) {
						
							// Pull in read length bases from reference to the beginning and end of the contig.
							prefix = c2r.getSequence(contigRead.getReferenceName(), 
									contigRead.getAlignmentStart()-readLength, readLength);
							suffix = c2r.getSequence(contigRead.getReferenceName(), contigRead.getAlignmentEnd()+1, readLength);
							
							bases = prefix.toUpperCase() + bases + suffix.toUpperCase();
						}
						
						Cigar cigar = new Cigar();
						if (contigRead.getCigarLength() == 1) {
							CigarElement elem = new CigarElement(first.getLength() + prefix.length() + suffix.length(), first.getOperator());
							cigar.add(elem);
						} else {
							CigarElement firstNew = new CigarElement(first.getLength() + prefix.length(), first.getOperator());
							CigarElement lastNew = new CigarElement(last.getLength() + suffix.length(), last.getOperator());
							
							cigar.add(firstNew);
							for (int i=1; i<contigRead.getCigarLength()-1; i++) {
								cigar.add(contigRead.getCigar().getCigarElement(i));
							}
							
							cigar.add(lastNew);
						}
						
						contigRead.setCigar(cigar);
						contigRead.setAlignmentStart(contigRead.getAlignmentStart()-prefix.length());

					} else {
						if (isDebug) {
							System.out.println("Not padding contig: " + contigRead.getReadName());
						}
					}
										
					//TODO: Safer delimiter?  This assumes no ~ in any read
					contigRead.setReadString("");
					String contigReadStr = contigRead.getSAMString();
					contigReadStr = contigReadStr.replace('\t','~');
					
					String contigName = contigRead.getReadName() + "_" + contigCount++ + "~" + contigReadStr; 
					
					writer.append(">" + contigName);
					writer.append(bases);
					writer.append("\n");
					hasCleanContigs = true;
				}
			}
		}
		contigReader.close();
		
		writer.close();
		
		return hasCleanContigs;
	}
	
	private String getPreprocessedBam(String tempDir) {
		return tempDir + "/" + "original_reads.bam";
	}
	
	private String getProprocessedFastq(String tempDir) {
		return tempDir + "/" + "original_reads.fastq.gz";
	}
	
	private String getTempReadFile(String tempDir) {
		if (isOutputIntermediateBam) {
			return getPreprocessedBam(tempDir);
		} else {
			return getProprocessedFastq(tempDir);
		}
	}

	SVReadCounter alignToSVContigs(String tempDir, String alignedToContigSam,
			String contigFasta, SAMFileWriter writer, SAMFileHeader header) throws IOException, InterruptedException {
		
		SVAlignerStdoutHandler stdoutHandler = new SVAlignerStdoutHandler(readLength, header);

		alignToContigs(tempDir, alignedToContigSam, contigFasta, writer, header, stdoutHandler);
		
		return stdoutHandler.getCounter();
	}
	
	void alignToContigs(String tempDir, String alignedToContigSam,
			String contigFasta, SAMFileWriter writer, SAMFileHeader header) throws IOException, InterruptedException {
		
		MutableBoolean isDone = new MutableBoolean();
		
		AdjustReadsQueueRunnable readQueueRunnable = new AdjustReadsQueueRunnable(threadManager, readAdjuster,
				writer, true, tempDir, header, isDone);
		
		AlignerStdoutHandler stdoutHandler = new AlignerStdoutHandler(readQueueRunnable);

		alignToContigs(tempDir, alignedToContigSam, contigFasta, writer, header, stdoutHandler);
	}
	
	void alignToContigs(String tempDir, String alignedToContigSam,
			String contigFasta, SAMFileWriter writer, SAMFileHeader header, StdoutHandler stdoutHandler) throws IOException, InterruptedException {
		
		String bam = getTempReadFile(tempDir);
		
		Aligner contigAligner = new Aligner(contigFasta, numThreads);
		
		// Align region fastq against assembled contigs
		contigAligner.shortAlign(bam, alignedToContigSam, stdoutHandler, isOutputIntermediateBam);
	}
	
	static class Pair<T, Y> {
		private T t;
		private Y y;
		public Pair(T t, Y y) {
			this.t = t;
			this.y = y;
		}
		
		public T getFirst() {
			return t;
		}
		
		public Y getSecond() {
			return y;
		}
	}
	
	static List<Feature> splitRegions(List<Feature> regions, 
			int maxRegionLength, int minRegionRemainder, int regionOverlap) {
		
		List<Feature> splitRegions = new ArrayList<Feature>();
		
		for (Feature region : regions) {
			if (region.getLength() <= maxRegionLength + minRegionRemainder) {
				splitRegions.add(region);
			} else {
				splitRegions.addAll(splitWithOverlap(region, maxRegionLength, minRegionRemainder, regionOverlap));
			}
		}
		
		return splitRegions;
	}
	
	/**
	 *  If any of the input list of features is greater than maxSize, split them into multiple features. 
	 */
	public static List<Feature> splitRegions(List<Feature> regions) {
		return splitRegions(regions, MAX_REGION_LENGTH, MIN_REGION_REMAINDER, REGION_OVERLAP);
	}
	
	public static List<Feature> splitWithOverlap(Feature region) {
		return splitWithOverlap(region, MAX_REGION_LENGTH, MIN_REGION_REMAINDER, REGION_OVERLAP);
	}
	
	static List<Feature> splitWithOverlap(Feature region, int maxRegionLength,
			int minRegionRemainder, int regionOverlap) {
		List<Feature> regions = new ArrayList<Feature>();
		
		long pos = region.getStart();
		long end = pos-1;
		
		while (end < region.getEnd()) {
			long start = pos;
			end = pos + maxRegionLength;
			long marker = end;
			
			// If we're at or near the end of the region, stop at region end.
			if (end > (region.getEnd() - minRegionRemainder)) {
				end = region.getEnd();
			}
			
			pos = marker - regionOverlap;
			
			regions.add(new Feature(region.getSeqname(), start, end));
		}
		
		return regions;
	}
	
	int[] getKmers(Feature region) {
		int[] kmerSizes = null;
		
		int kmerSize = region.getKmer();
		
		if (kmerSize > 0) {
			kmerSizes = toKmerArray(kmerSize, readLength);
		} else {
			kmerSizes = assemblerSettings.getKmerSize();
		}
		
		return kmerSizes;
	}
	
	int[] toKmerArray(int kmerSize, int readLength) {
		int[] kmerSizes = null;
		
		int maxKmerSize = this.readLength-15; 
		List<Integer> kmers = new ArrayList<Integer>();
		
		while (kmerSize < maxKmerSize) {
			kmers.add(kmerSize);
			kmerSize += 2;
		}
		
		kmerSizes = new int[kmers.size()];
		
		int i=0;
		for (int kmer : kmers) {
			kmerSizes[i++] = kmer;
		}

		return kmerSizes;
	}
		
	private NativeAssembler newAssembler(Feature region) {
		NativeAssembler assem = new NativeAssembler();

		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(assemblerSettings
				.getMaxPotentialContigs());

		assem.setMaxPathsFromRoot(100000);
		assem.setReadLength(readLength);
		//assem.setKmer(assemblerSettings.getKmerSize());
		assem.setKmer(getKmers(region));
		assem.setMinKmerFrequency(assemblerSettings.getMinNodeFrequncy());
		assem.setMinEdgeRatio(assemblerSettings.getMinEdgeRatio());
		assem.setMinBaseQuality(assemblerSettings.getMinBaseQuality());
		assem.setMinReadCandidateFraction(assemblerSettings.getMinReadCandidateFraction());
		assem.setMaxAverageDepth(assemblerSettings.getMaxAverageDepth());
		assem.setShouldSearchForSv(this.isPairedEnd && assemblerSettings.searchForStructuralVariation());
		assem.setAverageDepthCeiling(assemblerSettings.getAverageDepthCeiling());
		assem.setDebug(assemblerSettings.isDebug());

		return assem;
	}
	
	/*
	private NativeAssembler newUnalignedAssembler(int mnfMultiplier) {
		//Assembler assem = new JavaAssembler();
		NativeAssembler assem = new NativeAssembler();

		assem.setMaxContigs(MAX_POTENTIAL_UNALIGNED_CONTIGS);
		assem.setTruncateOutputOnRepeat(false);
		assem.setMaxPathsFromRoot(5000);
		assem.setReadLength(readLength);
		// Could be smaller for higher sensitivity here?
		int[] unalignedKmer = new int[1];
		unalignedKmer[0] = assemblerSettings.getKmerSize()[0];
		assem.setKmer(unalignedKmer);
		assem.setMinKmerFrequency(assemblerSettings.getMinUnalignedNodeFrequency());
		assem.setMinBaseQuality(assemblerSettings.getMinBaseQuality());

		return assem;
	}*/

	private void init() throws IOException {
		
		String javaVersion = System.getProperty("java.version");
		
		
		File workingDir = new File(tempDir);
		if (workingDir.exists()) {
			if (!workingDir.delete()) {
				throw new IllegalStateException("Unable to delete: " + tempDir);
			}
		}

		if (!workingDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir);
		}
		
		File unalignedTempDir = new File(tempDir + "/unaligned");
		
		if (!unalignedTempDir.mkdir()) {
			throw new IllegalStateException("Unable to create: " + tempDir + "/unaligned");
		}
		
		new NativeLibraryLoader().load(tempDir);
		
		threadManager = new ThreadManager(numThreads);
		
		if (this.localRepeatFile != null) {
			localRepeatWriter = new BufferedWriter(new FileWriter(localRepeatFile, false));
		}
	}
	
	private void mkdir(String dir) {
		File directory = new File(dir);
		if (!directory.mkdir()) {
			throw new IllegalStateException("Unable to create: " + dir);
		}
	}

	public void setReference(String reference) {
		this.reference = reference;
	}
	
	public void setBwaIndex(String bwaIndex) {
		this.bwaIndex = bwaIndex;
	}

	public void setTempDir(String temp) {
		this.tempDir = temp;
	}

	public void setAssemblerSettings(AssemblerSettings settings) {
		this.assemblerSettings = settings;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public void setMinContigMapq(int minContigMapq) {
		this.minContigMapq = minContigMapq;
	}
	
	public void setShouldReprocessUnaligned(boolean shouldReprocessUnaligned) {
		this.shouldReprocessUnaligned = shouldReprocessUnaligned;
	}
	
	public void setMaxUnalignedReads(int maxUnalignedReads) {
		this.maxUnalignedReads = maxUnalignedReads;
	}
	
	public CompareToReference2 getC2r() {
		return this.c2r;
	}
	
	public int getMinMappingQuality() {
		return this.minMappingQuality;
	}
	
	public int getMaxInsertLength() {
		return this.maxInsertLength;
	}
	
	public int getMinInsertLength() {
		return this.minInsertLength;
	}
	
	public void setMaxInsertLength(int maxInsertLen) {
		this.maxInsertLength = maxInsertLen;
	}
	
	public void setMinInsertLength(int minInsertLen) {
		this.minInsertLength = minInsertLen;
	}
	
	boolean isFiltered(SAMRecord read) {
		return SAMRecordUtils.isFiltered(isPairedEnd, read);
	}

	public static void run(String[] args) throws Exception {
		
		System.out.println("Starting 0.94 ...");
		
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(args);

		if (options.isValid()) {

			AssemblerSettings assemblerSettings = new AssemblerSettings();

			assemblerSettings.setKmerSize(options.getKmerSizes());
			assemblerSettings.setMinContigLength(options.getMinContigLength());
			assemblerSettings.setMinNodeFrequncy(options.getMinNodeFrequency());
			assemblerSettings.setMaxPotentialContigs(options
					.getMaxPotentialContigs());
			assemblerSettings.setMinUnalignedNodeFrequency(options.getMinUnalignedNodeFrequency());
			assemblerSettings.setMinBaseQuality(options.getMinBaseQuality());
			assemblerSettings.setMinReadCandidateFraction(options.getMinReadCandidateFraction());
			assemblerSettings.setMaxAverageDepth(options.getMaxAverageRegionDepth());
			assemblerSettings.setSearchForStructuralVariation(options.shouldSearchForStructuralVariation());
			assemblerSettings.setAverageDepthCeiling(options.getAverageDepthCeiling());
			assemblerSettings.setMinEdgeRatio(options.getMinEdgeRatio());
			assemblerSettings.setDebug(options.isDebug());

			ReAligner realigner = new ReAligner();
			realigner.setReference(options.getReference());
			realigner.setBwaIndex(options.getBwaIndex());
			realigner.setRegionsBed(options.getTargetRegionFile());
			realigner.setTempDir(options.getWorkingDir());
			realigner.setAssemblerSettings(assemblerSettings);
			realigner.setNumThreads(options.getNumThreads());
			realigner.setMinContigMapq(options.getMinContigMapq());
			realigner.setShouldReprocessUnaligned(!options.isSkipUnalignedAssembly());
			realigner.setMaxUnalignedReads(options.getMaxUnalignedReads());
			realigner.isPairedEnd = options.isPairedEnd();
			realigner.rnaSam = options.getRnaSam();
			realigner.rnaOutputSam = options.getRnaSamOutput();
			realigner.structuralVariantFile = options.getStructuralVariantFile();
			realigner.localRepeatFile = options.getLocalRepeatFile();
			realigner.minMappingQuality = options.getMinimumMappingQuality();
			realigner.hasPresetKmers = options.hasPresetKmers();
			realigner.isOutputIntermediateBam = options.useIntermediateBam();
			realigner.isDebug = options.isDebug();

			long s = System.currentTimeMillis();
			
			realigner.reAlign(options.getInputFiles(), options.getOutputFiles());

			long e = System.currentTimeMillis();

			System.out.println("Elapsed seconds: " + (e - s) / 1000);
		} else {
			System.exit(-1);
		}
	}
	
	public static void main(String[] args) throws Exception {
//		String inp = "--in /home/lmose/dev/ayc/opt/mem/test_tumor.bam --kmer 43 --mc-mapq 25 --mcl 101 --mcr -1.0 --mnf 2 --umnf 2 --mpc 50000 --out /home/lmose/dev/ayc/opt/mem/test_tumor.abra.bam --ref /home/lmose/reference/test/test.fa --targets /home/lmose/dev/ayc/opt/mem/test.gtf --threads 2 --working /home/lmose/dev/ayc/opt/mem/work1 --mur 50000000 --no-unalign --mbq 20 --rcf .02";
		String inp = "--in /home/lmose/dev/ayc/opt/mem/test_tumor.bam --kmer 43 --out /home/lmose/dev/ayc/opt/mem/test_tumor.abra3.bam --ref /home/lmose/reference/test/test.fa --targets /home/lmose/dev/ayc/opt/mem/test2.bed --threads 2 --working /home/lmose/dev/ayc/opt/mem/work3";
		run(inp.split("\\s+"));
	}
}

