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
import java.util.Random;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.picard.sam.SortSam;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * ABRA's main entry point
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReAligner {

	private static final int DEFAULT_MAX_UNALIGNED_READS = 1000000;
	
//	public static final int MAX_REGION_LENGTH = 2000;
//	private static final int MIN_REGION_REMAINDER = 500;
//	private static final int REGION_OVERLAP = 500;
	
	public static final int MAX_REGION_LENGTH = 400;
	private static final int MIN_REGION_REMAINDER = 200;
	private static final int REGION_OVERLAP = 200;

	private static final long RANDOM_SEED = 1;
	private static final int MAX_POTENTIAL_UNALIGNED_CONTIGS = 2000000;
	
	// Minimum sequence length recommended for use with bwa mem
	private static final int MIN_CONTIG_LENGTH = 70;
	
	private SAMFileHeader[] samHeaders;
	
	private List<Feature> regions;

	private String regionsGtf;

	private String tempDir;
	
	private String unalignedRegionSam;

	private String reference;
	
	private int minContigMapq;

	private AssemblerSettings assemblerSettings;
	
	private int numThreads;
	
	private int maxUnalignedReads = DEFAULT_MAX_UNALIGNED_READS;
	
	private boolean shouldReprocessUnaligned = true;
	
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
	
	private CompareToReference2 c2r;
	
	private ReadAdjuster readAdjuster;
	
	private ThreadManager threadManager;
	
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
		
		tempDirs = new String[inputSams.length];
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		writerFactory.setUseAsyncIo(true);
		
		writers = new SAMFileWriter[inputSams.length];
		
		for (int i=0; i<inputSams.length; i++) {
			// init temp dir
			String temp = tempDir + "/temp" + (i+1);
			mkdir(temp);
			tempDirs[i] = temp;

			// init BAM writer
			writers[i] = writerFactory.makeSAMOrBAMWriter(
					samHeaders[i], false, new File(outputFiles[i]));
		}
		
		// Start pre-processing reads on separate thread for each input file.
		// This happens in parallel with assembly, provided there are enough threads.
		for (int i=0; i<inputSams.length; i++) {
			preProcessReads(inputSams[i], tempDirs[i], writers[i]);
		}
		
		log("Iterating over regions");
		
		for (Feature region : regions) {
			log("Processing region: " + region.getDescriptor());
			spawnRegionThread(region, null);
		}
		
		log("Waiting for all threads to complete");
		threadManager.waitForAllThreadsToComplete();
		
		contigWriter.close();
		clock.stopAndPrint();
		
		String cleanContigsFasta = alignAndCleanContigs(contigFasta, tempDir, true);
				
		if (cleanContigsFasta != null) {		
			clock = new Clock("Align to contigs");
			clock.start();
			
			String[] alignedSams = alignReads(cleanContigsFasta, c2r);
			
			clock.stopAndPrint();

			clock = new Clock("Adjust reads");
			clock.start();
			log("Adjust reads");
			
			adjustReads(alignedSams, true, c2r);
			
			for (SAMFileWriter writer : this.writers) {
				writer.close();
			}
			
			clock.stopAndPrint();
			
			if (rnaSam != null) {
				processRna();
			}
			
		} else {
			log("WARNING!  No contigs assembled.  Just making a copy of input converting to/from SAM/BAM as appropriate.");
			for (int i=0; i<inputFiles.length; i++) {
				copySam(inputFiles[i], outputFiles[i]);	
			}
		}
		
		System.out.println("Done.");
	}
	
	private void preProcessReads(String inputSam, String tempDir, SAMFileWriter writer) throws InterruptedException {
		PreprocessReadsRunnable thread = new PreprocessReadsRunnable(threadManager, this,
				inputSam, this.getPreprocessedFastq(tempDir), c2r, writer);
		
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
		
		System.out.println("regions: " + regionsGtf);
		System.out.println("reference: " + reference);
		System.out.println("working dir: " + tempDir);
		System.out.println("num threads: " + numThreads);
		System.out.println("max unaligned reads: " + maxUnalignedReads);
		System.out.println(assemblerSettings.getDescription());
		System.out.println("rna: " + rnaSam);
		System.out.println("rna output: " + rnaOutputSam);
		System.out.println("paired end: " + isPairedEnd);
		
		System.out.println("Java version: " + System.getProperty("java.version"));
		
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
		Thread[] threads = new Thread[inputSams.length];
		
		for (int i=0; i<inputSams.length; i++) {
			alignedToContigsSams[i] = tempDirs[i] + "/" + "align_to_contig.sam";
			
			AlignReadsRunnable runnable = new AlignReadsRunnable(threadManager, this, tempDirs[i], inputSams[i], cleanContigsFasta,
					c2r, writers[i], alignedToContigsSams[i]);
			threads[i] = threadManager.spawnThread(runnable);
		}
		
		for (Thread thread : threads) {
			thread.join();
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

	private void adjustReads(String[] sortedAlignedToContig, 
			boolean isTightAlignment, CompareToReference2 c2r) throws InterruptedException, IOException {
		
		Thread[] threads = new Thread[inputSams.length];
		for (int i=0; i<inputSams.length; i++) {
			AdjustReadsRunnable runnable = new AdjustReadsRunnable(threadManager, readAdjuster, 
					sortedAlignedToContig[i], writers[i], isTightAlignment, tempDirs[i], samHeaders[i]);
			threads[i] = threadManager.spawnThread(runnable);
		}
		
		for (Thread thread : threads) {
			thread.join();
		}
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
			
//			System.out.println("chr: " + regionChromosome + " start: " + regionStart + "stop: " + regionStop);
			
//			String regionChromosome = fields[0];
//			int regionStart = Integer.parseInt(fields[1]) - 1000;
//			int regionStop = Integer.parseInt(fields[2]) + 1000;
			
			if ((contig.getReferenceName().equals(regionChromosome)) &&
				(contig.getAlignmentStart() >= regionStart) &&
				(contig.getAlignmentEnd() <= regionStop)) {
			
//				// Remove XP tags and other attributes, the semi-colons interfere with downstream processing
//				contig.clearAttributes();
//				contig.setAttribute("XP", null);
				
				outputReadsBam.addAlignment(contig);
			} else {
				System.out.println("Discarding: " + contig);
			}
		}
		
		outputReadsBam.close();
		reader.close();
	}
	
	private String alignAndCleanContigs(String contigFasta, String tempDir, boolean isTightAlignment) throws InterruptedException, IOException {
		log("Aligning contigs");
		Aligner aligner = new Aligner(reference, numThreads);
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
		cc.combine(contigsSam, contigsWithChim, isTightAlignment ? slack : 0);
		
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
			CompareToReference2 c2r, SAMFileWriter finalOutputSam, String alignedToContigSam) throws InterruptedException, IOException {
		log("Aligning original reads to contigs");
		alignToContigs(tempDir, inputSam, alignedToContigSam, cleanContigsFasta, c2r, finalOutputSam);
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
	
	private void downsampleSam(String sam, String downsampledSam, double keepProbability) {
		System.out.println("keepProbability: " + keepProbability);
		
		SAMFileReader reader = new SAMFileReader(new File(sam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter downsampleOutput = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeaders[0], true, new File(downsampledSam));

		Random random = new Random(RANDOM_SEED);
		int downsampleCount = 0;
		
		for (SAMRecord read : reader) {
			if (random.nextDouble() < keepProbability) {
				downsampleOutput.addAlignment(read);
				downsampleCount += 1;
			}
		}
		
		downsampleOutput.close();
		reader.close();
		
		System.out.println("Downsampled to: " + downsampleCount);
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
	
	/*
	private String getUnalignedReads(String unalignedBam) throws InterruptedException, IOException {
		
		int numUnalignedReads = 0;
		
		SAMFileWriter unalignedReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(unalignedBam));

		SAMFileReader reader = new SAMFileReader(new File(inputSam1));
		reader.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord read : reader) {
			if (shouldIncludeInUnalignedPile(read)) {
				unalignedReadsBam.addAlignment(read);
				numUnalignedReads += 1;
			}
		}
		
		reader.close();

		if (inputSam2 != null) {
			reader = new SAMFileReader(new File(inputSam2));
			reader.setValidationStringency(ValidationStringency.SILENT);
	
			for (SAMRecord read : reader) {
				if (shouldIncludeInUnalignedPile(read)) {
					unalignedReadsBam.addAlignment(read);
					numUnalignedReads += 1;
				}
			}
			
			reader.close();
		}
		
		if (inputSam3 != null) {
			reader = new SAMFileReader(new File(inputSam3));
			reader.setValidationStringency(ValidationStringency.SILENT);
	
			for (SAMRecord read : reader) {
				if (shouldIncludeInUnalignedPile(read)) {
					unalignedReadsBam.addAlignment(read);
					numUnalignedReads += 1;
				}
			}
			
			reader.close();
		}
		
		unalignedReadsBam.close();
		
		System.out.println("Number of unaligned reads: " + numUnalignedReads);
		
		if (numUnalignedReads > maxUnalignedReads) {
			double keepProbability = (double)  maxUnalignedReads / (double) numUnalignedReads;
			String downsampledSam = unalignedBam + ".downsampled.bam";
			downsampleSam(unalignedBam, downsampledSam, keepProbability);
			unalignedBam = downsampledSam;
		}
		
		return unalignedBam;		
	}
	*/
	
	private synchronized void appendContigs(String contigs) throws IOException {
		contigWriter.write(contigs);
	}
	
	public void processRegion(Feature region) throws Exception {
		
		try {
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			
			List<String> bams = new ArrayList<String>(Arrays.asList(this.inputSams));
			if (shouldReprocessUnaligned) {
				bams.add(unalignedRegionSam);
			}
			
			// Assemble contigs
			Assembler assem = newAssembler();
			String contigs = assem.assembleContigs(bams, contigsFasta, tempDir, region, region.getDescriptor(), true, this, c2r);
			if (!contigs.equals("<ERROR>") && !contigs.equals("<REPEAT>") && !contigs.isEmpty()) {
				appendContigs(contigs);
			}
		}
		catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}
		
	private void loadRegions() throws IOException {
		RegionLoader loader = new RegionLoader();
		regions = loader.load(regionsGtf);
		
//		RegionTracker regionTracker = new RegionTracker(regions, null);
//		regions = regionTracker.identifyTargetRegions(inputBams, this.assemblerSettings.getMinBaseQuality(), readLength, c2r);
		
		regions = collapseRegions(regions, readLength);
		regions = splitRegions(regions);
		
		System.out.println("Num regions: " + regions.size());
		for (Feature region : regions) {
			System.out.println(region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd());
		}
	}

	public void setRegionsGtf(String gtfFile) {
		this.regionsGtf = gtfFile;
	}

	private void getSamHeaderAndReadLength() {
		
		log("Identifying header and determining read length");
		
		this.samHeaders = new SAMFileHeader[this.inputSams.length];
		
		SAMFileReader reader = new SAMFileReader(new File(inputSams[0]));
		try {
			reader.setValidationStringency(ValidationStringency.SILENT);
	
			//ASSUMPTION!: All samples have same read length & insert len!
			samHeaders[0] = reader.getFileHeader();
			samHeaders[0].setSortOrder(SAMFileHeader.SortOrder.unsorted);
			
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
			
			System.out.println("Min insert length: " + minInsertLength);
			System.out.println("Max insert length: " + maxInsertLength);
			
		} finally {
			reader.close();
		}
		
		for (int i=1; i<this.inputSams.length; i++) {
			SAMFileReader reader2 = new SAMFileReader(new File(inputSams[i]));
			reader2.setValidationStringency(ValidationStringency.SILENT);
			samHeaders[i] = reader2.getFileHeader();
			samHeaders[i].setSortOrder(SAMFileHeader.SortOrder.unsorted);
			reader2.close();			
		}
				
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
	
	void sam2Fastq(String bam, String fastq, CompareToReference2 c2r, SAMFileWriter finalOutputSam) throws IOException {
		log("Preprocessing: " + bam);
		Sam2Fastq sam2Fastq = new Sam2Fastq();
		sam2Fastq.convert(bam, fastq, c2r, finalOutputSam, isPairedEnd, regions);
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

						// Pull in read length bases from reference to the beginning and end of the contig.
						String prefix = c2r.getSequence(contigRead.getReferenceName(), 
								contigRead.getAlignmentStart()-readLength, readLength);
						String suffix = c2r.getSequence(contigRead.getReferenceName(), contigRead.getAlignmentEnd()+1, readLength);
						
						bases = prefix.toUpperCase() + bases + suffix.toUpperCase();
						
						Cigar cigar = new Cigar();
						if (contigRead.getCigarLength() == 1) {
							// TODO: Pad here?
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
						System.out.println("Not padding contig: " + contigRead.getReadName());
					}
					
					// Aligned contigs are already expressed in forward strand context
	//				if (contigRead.getReadNegativeStrandFlag()) {
	//					bases = reverseComplementor.reverseComplement(bases);
	//				}
					
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
	
	private String getPreprocessedFastq(String tempDir) {
		return tempDir + "/" + "original_reads.fastq.gz";
	}
	
	private void alignToContigs(String tempDir, String inputSam, String alignedToContigSam,
			String contigFasta, CompareToReference2 c2r, SAMFileWriter finalOutputSam) throws IOException, InterruptedException {
		
		// Convert original bam to fastq
//		String fastq = tempDir + "/" + "original_reads.fastq.gz";
		String fastq = getPreprocessedFastq(tempDir);
		
		//TODO: Manage threads more intelligently based upon number of samples being processed.
		Aligner contigAligner = new Aligner(contigFasta, numThreads/2);
		
		// Align region fastq against assembled contigs
		contigAligner.shortAlign(fastq, alignedToContigSam);
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
	
	public static List<Feature> collapseRegions(List<Feature> regions, int readLength) {
		List<Feature> collapsedRegions = new ArrayList<Feature>();
		
		Feature currentRegion = null;
		
		for (Feature region : regions) {
			if (currentRegion != null) {
				if ((currentRegion.getSeqname().equals(region.getSeqname())) && 
					(currentRegion.getEnd() + (readLength) >= region.getStart())) {
					
					currentRegion.setEnd(region.getEnd());
				} else {
					collapsedRegions.add(currentRegion);
					currentRegion = region;
				}
			} else {
				currentRegion = region;
			}
		}
		
		if (currentRegion != null) {
			collapsedRegions.add(currentRegion);
		}
		
		System.out.println("Collapsed regions from " + regions.size() + " to " + collapsedRegions.size());
		
		return collapsedRegions;
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
			
//			if (end < region.getEnd()) {
//				end += regionOverlap;
//			}
			
			// If we're at or near the end of the region, stop at region end.
			if (end > (region.getEnd() - minRegionRemainder)) {
				end = region.getEnd();
			}
			
			pos = marker - regionOverlap;
			
			regions.add(new Feature(region.getSeqname(), start, end));
		}
		
		return regions;
	}
		
	private Assembler newAssembler() {
		NativeAssembler assem = new NativeAssembler();

		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(assemblerSettings
				.getMaxPotentialContigs());

		assem.setMaxPathsFromRoot(100000);
		assem.setReadLength(readLength);
		assem.setKmer(assemblerSettings.getKmerSize());
		assem.setMinKmerFrequency(assemblerSettings.getMinNodeFrequncy());
		assem.setMinBaseQuality(assemblerSettings.getMinBaseQuality());
		assem.setMinReadCandidateFraction(assemblerSettings.getMinReadCandidateFraction());

		return assem;
	}
	
	private Assembler newUnalignedAssembler(int mnfMultiplier) {
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
	}

	private void init() throws IOException {
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
		
		System.out.println("Starting 0.74 ...");
		
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

			ReAligner realigner = new ReAligner();
			realigner.setReference(options.getReference());
			realigner.setRegionsGtf(options.getTargetRegionFile());
			realigner.setTempDir(options.getWorkingDir());
			realigner.setAssemblerSettings(assemblerSettings);
			realigner.setNumThreads(options.getNumThreads());
			realigner.setMinContigMapq(options.getMinContigMapq());
			realigner.setShouldReprocessUnaligned(!options.isSkipUnalignedAssembly());
			realigner.setMaxUnalignedReads(options.getMaxUnalignedReads());
			realigner.isPairedEnd = options.isPairedEnd();
			realigner.rnaSam = options.getRnaSam();
			realigner.rnaOutputSam = options.getRnaSamOutput();

			long s = System.currentTimeMillis();
			
			realigner.reAlign(options.getInputFiles(), options.getOutputFiles());

			long e = System.currentTimeMillis();

			System.out.println("Elapsed seconds: " + (e - s) / 1000);
		}
	}
	
	public static void main(String[] args) throws Exception {
		String inp = "--in /home/lmose/dev/ayc/opt/mem/test_tumor.bam --kmer 43 --mc-mapq 25 --mcl 101 --mcr -1.0 --mnf 2 --umnf 2 --mpc 50000 --out /home/lmose/dev/ayc/opt/mem/test_tumor.abra.bam --ref /home/lmose/reference/test/test.fa --targets /home/lmose/dev/ayc/opt/mem/test.gtf --threads 2 --working /home/lmose/dev/ayc/opt/mem/work1 --mur 50000000 --no-unalign --mbq 20 --rcf .02";
		run(inp.split("\\s+"));
	}
	/*
	public static void main(String[] args) throws Exception {
		ReAligner realigner = new ReAligner();
//		String originalReadsSam = args[0];
//		String alignedToContigSam = args[1];
//		String unalignedSam = args[2];
//		String outputFilename = args[3];
		
//		String originalReadsSam = "/home/lmose/dev/ayc/sim/s43/orig_atc.bam";
//		String alignedToContigSam = "/home/lmose/dev/ayc/sim/s43/atc.bam";
//		String unalignedSam = args[2];
//		String outputFilename = "/home/lmose/dev/ayc/sim/s43/atc_out.bam";
		
		SAMFileReader reader = new SAMFileReader(new File("/home/lmose/dev/abra/missing_68i/header.sam"));
	
		realigner.samHeader = reader.getFileHeader();
		
		reader.close();
		
//		realigner.readLength = 76;
//		realigner.isPairedEnd = true;
		realigner.minInsertLength = 100;
		realigner.maxInsertLength = 500;
//		realigner.regionsGtf = "/home/lmose/dev/ayc/p3/wxs.gtf";
//		realigner.loadRegions();
		
//		realigner.getSamHeader(originalReadsSam);
		
		realigner.isPairedEnd = true;
				
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("/home/lmose/reference/chr21/21.fa");
		realigner.c2r = c2r;
		
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				realigner.samHeader, false, new File("/home/lmose/dev/abra/missing_68i/output.bam"));

//		ReadAdjuster readAdjuster = new ReadAdjuster(realigner, realigner.maxMapq);
//		readAdjuster.adjustReads("/home/lmose/dev/abra/missing_68i/a2c.sam", writer, true, c2r, "/home/lmose/dev/abra/missing_68i/temp", realigner.samHeader);
		
		writer.close();
	}
*/

	/*
	public static void main(String[] args) throws Exception {
		System.out.println("Adjusting 2...");
		ReAligner ra = new ReAligner();
		ra.isPairedEnd = false;
		
//		String headerSourceBam = args[0];
//		String alignedToContigBam = args[1];
//		String outputBam = args[2];
//		String reference = args[3];
//		String tempDir = args[4];
		
//		String headerSourceBam = "/home/lmose/dev/ayc/paired/header.bam";
//		String alignedToContigBam = "/home/lmose/dev/ayc/paired/a2c.bam";
//		String outputBam = "/home/lmose/dev/ayc/paired/output.bam";
//		String tempDir = "/home/lmose/dev/ayc/paired/adjust_temp";
//		String reference = args[3];
		
		String headerSourceBam = "/home/lmose/dev/ayc/unaligned/chr1.bam";
		String alignedToContigBam = "/home/lmose/dev/ayc/unaligned/align_to_contig.bam";
		String outputBam = "/home/lmose/dev/ayc/unaligned/output.bam";
		String tempDir = "/home/lmose/dev/ayc/unaligned/adjust_temp";
		String reference = "/home/lmose/reference/chr1b/chr1.fa";

		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
//		writerFactory.setUseAsyncIo(true);
		
		ra.getSamHeaderAndReadLength(headerSourceBam);
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				ra.samHeader, false, new File(outputBam));
				
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(reference);
		
//		CompareToReference2 c2r = null;
		
		ra.adjustReads(alignedToContigBam, writer, false, c2r, tempDir);
				
		writer.close();
		
		System.out.println("Done");
	}
	*/
	
	
	/*
	public void testA2c() throws Exception {
		SAMFileReader reader = new SAMFileReader(new File("/home/lmose/dev/ayc/40c/align_to_contig.bam"));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		int count = 0;
		
		SamStringReader samStringReader = new SamStringReader(reader.getFileHeader());
		for (SAMRecord read : reader) {
			
			try {
				int numBestHits = getIntAttribute(read, "X0");
				int numSubOptimalHits = getIntAttribute(read, "X1");
				
				int totalHits = numBestHits + numSubOptimalHits;
				
				//TODO: If too many best hits, what to do?
				
	//			spikeLog("total hits: " + totalHits, origRead);
				
				if ((totalHits > 1) && (totalHits < 1000)) {
					String alternateHitsStr = (String) read.getAttribute("XA");
						
						String[] alternates = alternateHitsStr.split(";");
						
						for (int i=0; i<alternates.length-1; i++) {
							String[] altInfo = alternates[i].split(",");
							String altContigReadStr = altInfo[0];
							char strand = altInfo[1].charAt(0);
							int position = Integer.parseInt(altInfo[1].substring(1));
							String cigar = altInfo[2];
							int mismatches = Integer.parseInt(altInfo[3]);
							
							altContigReadStr = altContigReadStr.substring(altContigReadStr.indexOf('~')+1);
							altContigReadStr = altContigReadStr.replace('~', '\t');
							SAMRecord contigRead = samStringReader.getRead(altContigReadStr);
						}
	
				}
				
				System.out.println("count: " + count);
				if ((count % 100000) == 0) {
					System.out.println("count: " + count);
				}
				count++;
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println(read.getSAMString());
				throw e;
			}
		}
		
		reader.close();
	}
	
	public static void main(String[] args) throws Exception {
		new ReAligner().testA2c();
	}
	*/
	
	/*
	public static void main(String[] args) throws Exception {
		System.out.println("0.2");
		ReAligner realigner = new ReAligner();

		long s = System.currentTimeMillis();

		// sim204 test
//		String input = "/home/lmose/dev/ayc/sim/s204/chr2.bam";
//		String input2 = "/home/lmose/dev/ayc/sim/s204/chr1.bam";
//		String output = "/home/lmose/dev/ayc/sim/s204/normal.abra.bam";
//		String output2 = "/home/lmose/dev/ayc/sim/s204/tumor.abra.bam";
//		String reference = "/home/lmose/reference/chr1b/chr1.fa";
////		String regions = "/home/lmose/dev/abra_wxs/4/4.gtf";
//		String regions = "/home/lmose/dev/ayc/sim/s204/204.gtf";
//		String tempDir = "/home/lmose/dev/ayc/sim/s204/working";

		
		// mem test
//		String input = "/home/lmose/dev/ayc/sim/s204/chr2.bam";
//		String input2 = "/home/lmose/dev/ayc/sim/s204/tchr1.bam";
//		String output = "/home/lmose/dev/ayc/sim/s204/normal.abra2.bam";
//		String output2 = "/home/lmose/dev/ayc/sim/s204/tumor.abra2.bam";
//		String reference = "/home/lmose/reference/chr1b/chr1.fa";
//		String regions = "/home/lmose/dev/ayc/regions/clinseq5/uncseq5_chr1.gtf";
//		String tempDir = "/home/lmose/dev/ayc/sim/s204/working2";
//
		
		String input = "/home/lmose/dev/ayc/sim/s204/chr2.bam";
		String input2 = "/home/lmose/dev/ayc/sim/s204/t2.bam";
		String output = "/home/lmose/dev/ayc/sim/s204/normal.abra3.bam";
		String output2 = "/home/lmose/dev/ayc/sim/s204/tumor.abra3.bam";
		String reference = "/home/lmose/reference/chr1b/chr1.fa";
		String regions = "/home/lmose/dev/ayc/sim/s204/204.gtf";
		String tempDir = "/home/lmose/dev/ayc/sim/s204/working3";



		AssemblerSettings settings = new AssemblerSettings();
		settings.setKmerSize(63);
		settings.setMinContigLength(100);
		settings.setMinEdgeFrequency(2);
		settings.setMinNodeFrequncy(2);
		settings.setMaxPotentialContigs(30000);
		settings.setMinContigRatio(-1.0);

		realigner.setAssemblerSettings(settings);
		
		realigner.setMinContigMapq(1);
		realigner.setReference(reference);
		realigner.setRegionsGtf(regions);
		realigner.setTempDir(tempDir);
//		realigner.setNumThreads(4);
		realigner.setNumThreads(2);

		realigner.reAlign(input, input2, output, output2);

		long e = System.currentTimeMillis();

		System.out.println("Elapsed seconds: " + (e - s) / 1000);
		
//		Thread.sleep(600000);
	}
	*/
}

