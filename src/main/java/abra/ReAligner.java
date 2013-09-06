/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

//import abra.Assembler.TooManyPotentialContigsException;

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
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class ReAligner {

	private static final int DEFAULT_MAX_UNALIGNED_READS = 1000000;
	public static final int MAX_REGION_LENGTH = 2000;
//	public static final int MAX_REGION_LENGTH = 25000;
	private static final int MIN_REGION_REMAINDER = 500;
//	private static final int MIN_REGION_REMAINDER = 300;
//	private static final int REGION_OVERLAP = 200; 
	private static final int REGION_OVERLAP = 500;
	private static final long RANDOM_SEED = 1;
	private static final int MAX_POTENTIAL_UNALIGNED_CONTIGS = 2000000;
	
	private int missingXATag = 0;
	
	private SAMFileHeader samHeader;

	private long startMillis;

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
	
	private List<ReAlignerRunnable> threads = new ArrayList<ReAlignerRunnable>();
	
	private ReverseComplementor reverseComplementor = new ReverseComplementor();
	
	private String inputSam1;
	private String inputSam2;
	private String inputSam3;
	
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
	
	public void reAlign(String inputSam, String inputSam2, String inputSam3, String outputSam, String outputSam2, String outputSam3) throws Exception {

		System.out.println("input: " + inputSam);
		System.out.println("input2: " + inputSam2);
		System.out.println("input3: " + inputSam3);
		System.out.println("output: " + outputSam);
		System.out.println("output2: " + outputSam2);
		System.out.println("output3: " + outputSam3);
		System.out.println("regions: " + regionsGtf);
		System.out.println("reference: " + reference);
		System.out.println("working dir: " + tempDir);
		System.out.println("num threads: " + numThreads);
		System.out.println("max unaligned reads: " + maxUnalignedReads);
		System.out.println(assemblerSettings.getDescription());
		System.out.println("rna: " + rnaSam);
		System.out.println("rna output: " + rnaOutputSam);
		
		System.out.println("Java version: " + System.getProperty("java.version"));

		startMillis = System.currentTimeMillis();

		this.inputSam1 = inputSam;
		this.inputSam2 = inputSam2;
		this.inputSam3 = inputSam3;
		
		init();
		
		System.out.println("c2r1: " + c2r);
		c2r = new CompareToReference2();
		c2r.init(this.reference);
		System.out.println("c2r2: " + c2r);

		log("Reading Input SAM Header and identifying read length");
		getSamHeaderAndReadLength(inputSam);
		
		log("Read length: " + readLength);
		
		log("Loading target regions");
		loadRegions();
		
		samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
		
		if (shouldReprocessUnaligned) {
			Clock clock = new Clock("Process unaligned");
			clock.start();
			log("Assembling unaligned reads");
			
			String unalignedSam = tempDir + "/unaligned.bam";
			unalignedSam = getUnalignedReads(unalignedSam);
			
//			String unalignedSam = tempDir + "/" + "unaligned_to_contig.bam";
			
			String unalignedDir = tempDir + "/unaligned";
			String unalignedContigFasta = unalignedDir + "/unaligned_contigs.fasta";
			unalignedRegionSam = unalignedDir + "/unaligned_region.bam";
			String sortedUnalignedRegion = unalignedDir + "/sorted_unaligned_region.bam";
			
			Assembler assem = newUnalignedAssembler(1);
			List<String> unalignedSamList = new ArrayList<String>();
			unalignedSamList.add(unalignedSam);
			boolean hasContigs = assem.assembleContigs(unalignedSamList, unalignedContigFasta, tempDir, null, "unaligned", false, this);

			// Make eligible for GC
			assem = null;
						
			if (hasContigs) {
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
					adjustReads(alignedToContigBam, unalignedWriter, false, null, unalignedDir);
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
		
		Clock clock = new Clock("Assembly");
		clock.start();
		
		String contigFasta = tempDir + "/" + "all_contigs.fasta";
		contigWriter = new BufferedWriter(new FileWriter(contigFasta, false));
		
		log("Iterating over regions");
		for (Feature region : regions) {
			log("Spawning thread for: " + region.getDescriptor());
//			spawnRegionThread(region, assemblyBam);
			spawnRegionThread(region, null);
		}
		
		log("Waiting for all threads to complete");
		waitForAllThreadsToComplete();
		
		contigWriter.close();
		clock.stopAndPrint();
		
//		clock = new Clock("Combine contigs");
//		clock.start();
//		log("Combining contigs");
//		
//		combineContigs(contigFasta);
//		clock.stopAndPrint();
		
		String cleanContigsFasta = alignAndCleanContigs(contigFasta, tempDir, true);
				
		if (cleanContigsFasta != null) {
			String tempDir1 = tempDir + "/temp1";
			String tempDir2 = tempDir + "/temp2";
			String tempDir3 = tempDir + "/temp3";
			mkdir(tempDir1);
			mkdir(tempDir2);
			mkdir(tempDir3);
			
			SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
//			writerFactory.setUseAsyncIo(true);
			
			SAMFileWriter writer1 = writerFactory.makeSAMOrBAMWriter(
					samHeader, false, new File(outputSam));
			
			SAMFileWriter writer2 = null;
			
			SAMFileWriter writer3 = null;
			
			if (inputSam2 != null) {
				writer2 = writerFactory.makeSAMOrBAMWriter(
					samHeader, false, new File(outputSam2));
			}
			
			if (inputSam3 != null) {
				writer3 = writerFactory.makeSAMOrBAMWriter(
					samHeader, false, new File(outputSam3));
			}
			
			clock = new Clock("Sam2Fastq and Align");
			clock.start();
//			String alignedToContigSam1 = alignReads(tempDir1, inputSam, cleanContigsFasta, c2r, writer1);
//			String alignedToContigSam2 = null;
//			if (inputSam2 != null) {
//				alignedToContigSam2 = alignReads(tempDir2, inputSam2, cleanContigsFasta, c2r, writer2);
//			}
			
			System.out.println("c2r3: " + c2r);
			
			String[] alignedSams = alignReads(tempDir1, inputSam, cleanContigsFasta, c2r, writer1, 
					tempDir2, inputSam2, writer2, tempDir3, inputSam3, writer3);
			
			String alignedToContigSam1 = alignedSams[0];
			String alignedToContigSam2 = alignedSams[1];
			String alignedToContigSam3 = alignedSams[2];
			
			clock.stopAndPrint();

			String alignedToContigBam1 = alignedToContigSam1;
			String alignedToContigBam2 = alignedToContigSam2;
			String alignedToContigBam3 = alignedToContigSam3;

			clock = new Clock("Adjust reads");
			clock.start();
			log("Adjust reads");
			// Output sorted by coordinate
			//samHeader.setSortOrder(SortOrder.coordinate);
			samHeader.setSortOrder(SortOrder.unsorted);
			
			adjustReads(alignedToContigBam1, writer1,
					alignedToContigBam2, writer2, true, c2r, tempDir1, tempDir2,
					alignedToContigBam3, writer3, tempDir3);
			
			writer1.close();
			if (writer2 != null) {
				writer2.close();
			}
			
			if (writer3 != null) {
				writer3.close();
			}
			
			clock.stopAndPrint();
			
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
		} else {
			log("WARNING!  No contigs assembled.  Just making a copy of input converting to/from SAM/BAM as appropriate.");
			copySam(inputSam, outputSam);
			if (inputSam2 != null) {
				copySam(inputSam2, outputSam2);
			}
			if (inputSam3 != null) {
				copySam(inputSam2, outputSam2);
			}
		}
		
		System.out.println("Done.");
	}
	
	private void copySam(String input, String output) {
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				samHeader, false, new File(output));
		
		for (SAMRecord read : reader) {
			writer.addAlignment(read);
		}
		
		reader.close();
		writer.close();
	}
	
	private String[] alignReads(String tempDir1, String inputSam, 
			String cleanContigsFasta, CompareToReference2 c2r, SAMFileWriter writer1, 
			String tempDir2, String inputSam2, SAMFileWriter writer2,
			String tempDir3, String inputSam3, SAMFileWriter writer3) throws IOException, InterruptedException {
		
		System.out.println("c2r4: " + c2r);
		
		// Build contig fasta index
		log("Indexing contigs");
		Aligner contigAligner = new Aligner(cleanContigsFasta, numThreads);
		contigAligner.index();
		log("Contig indexing done");
		
		
		String alignedToContigSam1 = tempDir1 + "/" + "align_to_contig.sam";
		AlignReadsRunnable alignReadsRunnable1 = new AlignReadsRunnable(this,
				tempDir1, inputSam, cleanContigsFasta, c2r, writer1, alignedToContigSam1);
		
		Thread thread1 = new Thread(alignReadsRunnable1);
		thread1.start();
		
//		alignReads(tempDir1, inputSam, cleanContigsFasta, c2r, writer1, alignedToContigSam1);
		
		String alignedToContigSam2 = null;
		
		Thread thread2 = null;
		
		if (inputSam2 != null) {
			alignedToContigSam2 = tempDir2 + "/" + "align_to_contig.sam";
			AlignReadsRunnable alignReadsRunnable2 = new AlignReadsRunnable(this,
					tempDir2, inputSam2, cleanContigsFasta, c2r, writer2, alignedToContigSam2);
			thread2 = new Thread(alignReadsRunnable2);
			thread2.start();
			
//			alignReads(tempDir2, inputSam2, cleanContigsFasta, c2r, writer2, alignedToContigSam2);
		}
		
		String alignedToContigSam3 = null;
		
		Thread thread3 = null;
		
		if (inputSam3 != null) {
			alignedToContigSam3 = tempDir3 + "/" + "align_to_contig.sam";
			AlignReadsRunnable alignReadsRunnable3 = new AlignReadsRunnable(this,
					tempDir3, inputSam3, cleanContigsFasta, c2r, writer3, alignedToContigSam3);
			thread3 = new Thread(alignReadsRunnable3);
			thread3.start();
		}
		
		thread1.join();
		
		if (inputSam2 != null) {
			thread2.join();
		}
		
		if (inputSam3 != null) {
			thread3.join();
		}
		
		return new String[] { alignedToContigSam1, alignedToContigSam2, alignedToContigSam3 };
	}
	
	/*
	public void findFusions(String inputSam, String outputSam) {
		this.inputSam1 = inputSam;
		
		init();
		
		String fusDir = tempDir + "/fusion";
		
		if (!new File(fusDir).mkdir()) {
			throw new RuntimeException("Failed to create: " + fusDir);
		}
		
		log("Reading Input SAM Header and identifying read length");
		getSamHeaderAndReadLength(inputSam);
		
		String candidateBam = fusDir + "/candidates.bam";
		
		SAMFileWriter candidateWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(candidateBam));
		
		SamReadPairReader reader = new SamReadPairReader(inputSam);
		
		for (ReadPair pair : reader) {
			SAMRecord read1 = pair.getRead1();
			SAMRecord read2 = pair.getRead2();
			
			boolean isCandidate = false;
			
			if (!read1.getReadUnmappedFlag() && !read2.getReadUnmappedFlag()) {
				if (!read1.getReferenceName().equals(read2.getReferenceName())) {
					isCandidate = true;
				} else if (Math.abs(read1.getAlignmentStart() - read2.getAlignmentStart()) > 2000) {
					isCandidate = true;
				}
			}
			
			if (isCandidate) {
				candidateWriter.addAlignment(read1);
				candidateWriter.addAlignment(read2);
			}
		}
		
		reader.close();
		candidateWriter.close();
	}
	*/
		
	void updateMismatchAndEditDistance(SAMRecord read, CompareToReference2 c2r, SAMRecord origRead) {
		if (read.getAttribute("YO") != null) {
			int numMismatches = c2r.numMismatches(read);				
			int numIndelBases = getNumIndelBases(read);
			read.setAttribute("XM", numMismatches);
			read.setAttribute("NM", numMismatches + numIndelBases);
			read.setMappingQuality(calcMappingQuality(read, origRead));
			
			if (numMismatches > (read.getReadLength()/10)) {
				System.out.println("HIGH_MISMATCH: [" + read.getSAMString() + "]");
			}
		}
	}
		
	//TODO: Add rhyme or reason to this
	private int calcMappingQuality(SAMRecord read, SAMRecord origRead) {
		int mapq = 0;
		
		// Need original read here because updated read has already had 0x04 flag unset.
		
		if ((origRead.getReadUnmappedFlag()) || (read.getMappingQuality() > 0)) {
			int contigQuality = (Integer) read.getAttribute("YQ");
			int quality = Math.min(contigQuality, this.maxMapq);
			int mismatchesToContig = (Integer) read.getAttribute("YM");
			quality -= mismatchesToContig * 5;
			mapq = Math.max(quality, 1);
		}
		
		return mapq;
	}
	
	public static int getNumIndels(SAMRecord read) {
		int numIndels = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I)) {
				numIndels += 1;
			}
		}
		
		return numIndels;
	}

	public static int getNumIndelBases(SAMRecord read) {
		int numIndelBases = 0;
		
		for (CigarElement element : read.getCigar().getCigarElements()) {
			if ((element.getOperator() == CigarOperator.D) || (element.getOperator() == CigarOperator.I)) {
				numIndelBases += element.getLength();
			}
		}
		
		return numIndelBases;
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

	private void adjustReads(String sortedAlignedToContig1, SAMFileWriter outputSam1, 
			String sortedAlignedToContig2, SAMFileWriter outputSam2, boolean isTightAlignment,
			CompareToReference2 c2r, String tempDir1, String tempDir2,
			String sortedAlignedToContig3, SAMFileWriter outputSam3, String tempDir3) throws InterruptedException, IOException {
		
		if (this.numThreads > 1) {
			
			System.out.println("Adjusting reads in parallel");
			AdjustReadsRunnable runnable1 = new AdjustReadsRunnable(this, sortedAlignedToContig1, outputSam1, isTightAlignment, c2r, tempDir1);
			Thread thread1 = new Thread(runnable1);
			thread1.start();
			
			Thread thread2 = null;
			if (inputSam2 != null) {
				AdjustReadsRunnable runnable2 = new AdjustReadsRunnable(this, sortedAlignedToContig2, outputSam2, isTightAlignment, c2r, tempDir2);
				thread2 = new Thread(runnable2);
				thread2.start();
			}
			
			Thread thread3 = null;
			if (inputSam3 != null) {
				AdjustReadsRunnable runnable3 = new AdjustReadsRunnable(this, sortedAlignedToContig3, outputSam3, isTightAlignment, c2r, tempDir3);
				thread3 = new Thread(runnable3);
				thread3.start();
			}
			
			thread1.join();
			
			if (inputSam2 != null) {
				thread2.join();
			}
			
			if (inputSam3 != null) {
				thread3.join();
			}
			
		} else {
			System.out.println("Adjusting reads sequentially");
			adjustReads(sortedAlignedToContig1, outputSam1, isTightAlignment, c2r, tempDir1);
			if (inputSam2 != null) {
				adjustReads(sortedAlignedToContig2, outputSam2, isTightAlignment, c2r, tempDir2);
			}
			if (inputSam3 != null) {
				adjustReads(sortedAlignedToContig3, outputSam3, isTightAlignment, c2r, tempDir3);
			}
		}
	}
	
	private void discardMisalignedContigs(String inputSam, String outputSam) {
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter outputReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(outputSam));

		for (SAMRecord contig : reader) {
			String[] fields = contig.getReadName().split("_");
			String regionChromosome = fields[0];
			int regionStart = Integer.parseInt(fields[1]) - 1000;
			int regionStop = Integer.parseInt(fields[2]) + 1000;
			
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
		aligner.align(contigFasta, contigsSam);
		
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
		
//			log("Loading reference for chopper clopper.");
//			CompareToReference2 c2r = new CompareToReference2();
//			c2r.init(this.reference);
			
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
	
	private void sortBams(String in1, String in2, String out1, String out2, String sortOrder) throws InterruptedException {
		if (numThreads > 1) {
			System.out.println("Sorting BAMs in parallel");
			SortBamRunnable runnable1 = new SortBamRunnable(this, in1, out1, sortOrder);
			Thread thread1 = new Thread(runnable1);
			
			SortBamRunnable runnable2 = new SortBamRunnable(this, in2, out2, sortOrder);
			Thread thread2 = new Thread(runnable2);

			thread1.start();
			thread2.start();
			
			thread1.join();
			thread2.join();
		} else {
			System.out.println("Sorting bams sequentially");
			sortBam(in1, out1, sortOrder);
			sortBam(in2, out2, sortOrder);
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
			
	private void combineContigs(String contigFasta) throws IOException, InterruptedException {

		System.out.println("Combining contigs...");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(contigFasta, false));
		
		for (Feature region : regions) {
			String regionContigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			BufferedReader reader = new BufferedReader(new FileReader(regionContigsFasta));
			String line = reader.readLine();
			while (line != null) {
				writer.write(line);
				writer.write('\n');
				line = reader.readLine();
			}
			reader.close();
		}
		
		writer.close();
	}
	
	public synchronized void addThread(ReAlignerRunnable thread) {
		threads.add(thread);
	}
	
	public synchronized void removeThread(ReAlignerRunnable thread) {
		threads.remove(thread);
	}
	
	private synchronized int activeThreads() {
		return threads.size();
	}
	
	private void waitForAvailableThread() throws InterruptedException {
		while (activeThreads() == numThreads) {
			Thread.sleep(500);
		}
	}
	
	private void waitForAllThreadsToComplete() throws InterruptedException, IOException {
		long start = System.currentTimeMillis();
		while (activeThreads() > 0) {
			long elapsedSecs = (System.currentTimeMillis() - start) / 1000;
			if ((elapsedSecs % 60) == 0) {
				log("Waiting on " + threads.size() + " threads.");
				
				logOSMemory();
				
			}
			Thread.sleep(500);
		}
	}
	
	private void spawnRegionThread(Feature region, String inputSam) throws InterruptedException {
		waitForAvailableThread();
//		ReAlignerRunnable thread = new ReAlignerRunnable(this, region, inputSam);
		ReAlignerRunnable thread = new ReAlignerRunnable(this, region, null);
		addThread(thread);
		new Thread(thread).start();
	}
	
	private void downsampleSam(String sam, String downsampledSam, double keepProbability) {
		System.out.println("keepProbability: " + keepProbability);
		
		SAMFileReader reader = new SAMFileReader(new File(sam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriter downsampleOutput = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(downsampledSam));

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
			/*
			// Consider long soft clipping
			else if (read.getCigarLength() > 1) {
				CigarElement first = read.getCigar().getCigarElement(0);
				CigarElement last = read.getCigar().getCigarElement(read.getCigarLength()-1);
			
				int clipLength = 0;
				if (first.getOperator() == CigarOperator.S) {
					clipLength += first.getLength();
				}
				if (last.getOperator() == CigarOperator.S) {
					clipLength += last.getLength(); 
				}
				
				if (clipLength >= 10) {
					shouldInclude = true;
				}
			}
			*/
		}
		
		return shouldInclude;
	}
	
	private String getUnalignedReads(String unalignedBam) throws InterruptedException, IOException {
		
		int numUnalignedReads = 0;
		
		SAMFileWriter unalignedReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				samHeader, true, new File(unalignedBam));

		SAMFileReader reader = new SAMFileReader(new File(inputSam1));
		reader.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord read : reader) {
			//TODO: Remove hardcoded mapq.  Exclude from original region?
//			if (read.getReadUnmappedFlag() || read.getMappingQuality() < 10) {
//		if (read.getReadUnmappedFlag()) {
//			if ((read.getReadUnmappedFlag()) || (read.getCigarLength() >= 4)) {
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
				//TODO: Remove hardcoded mapq.  Exclude from original region?
	//			if (read.getReadUnmappedFlag() || read.getMappingQuality() < 10) {
	//			if (read.getReadUnmappedFlag()) {
	//			if ((read.getReadUnmappedFlag()) || (read.getCigarLength() >= 4)) {
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
				//TODO: Remove hardcoded mapq.  Exclude from original region?
	//			if (read.getReadUnmappedFlag() || read.getMappingQuality() < 10) {
	//			if (read.getReadUnmappedFlag()) {
	//			if ((read.getReadUnmappedFlag()) || (read.getCigarLength() >= 4)) {
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
	
	private synchronized void appendContigs(BufferedReader reader) throws IOException {
		long start = System.currentTimeMillis();
		
		String line = reader.readLine();
		while (line != null) {
			contigWriter.write(line);
			contigWriter.write('\n');
			line = reader.readLine();
		}
		
		long end = System.currentTimeMillis();
		
		System.out.println("Elapsed msecs in appendConting: " + (end-start));
	}
	
	public void processRegion(Feature region) throws Exception {
		
		try {
			String contigsFasta = tempDir + "/" + region.getDescriptor() + "_contigs.fasta";
			
			List<String> bams = new ArrayList<String>();
			bams.add(inputSam1);
			if (inputSam2 != null) {
				bams.add(inputSam2);
			}
			if (inputSam3 != null) {
				bams.add(inputSam3);
			}
			if (shouldReprocessUnaligned) {
				bams.add(unalignedRegionSam);
			}
			
			// Assemble contigs
			Assembler assem = newAssembler();
			boolean hasContigs = assem.assembleContigs(bams, contigsFasta, tempDir, region, region.getDescriptor(), true, this);
			
			// Append contigs to the global fasta file
			if (hasContigs) {
				BufferedReader reader = new BufferedReader(new FileReader(contigsFasta));
				appendContigs(reader);
				reader.close();
			}
			
			// Now delete the temporary fasta file.
			File localAssembledContigs = new File(contigsFasta);
			localAssembledContigs.delete();
		}
		catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}
	
	private void loadRegions() throws IOException {
		GtfLoader loader = new GtfLoader();
		regions = loader.load(regionsGtf);
		
		regions = collapseRegions(regions, readLength);
		
		regions = splitRegions(regions);
		
		for (Feature region : regions) {
			System.out.println(region.getSeqname() + "\t" + region.getStart() + "\t" + region.getEnd());
		}
	}

	public void setRegionsGtf(String gtfFile) {
		this.regionsGtf = gtfFile;
	}

	private void log(String message) {
		System.out.println(new Date() + " : " + message);
	}

	private void getSamHeaderAndReadLength(String inputSam) {
		
		log("Identifying header and determining read length");
		SAMFileReader reader = new SAMFileReader(new File(inputSam));
		try {
			reader.setValidationStringency(ValidationStringency.SILENT);
	
			samHeader = reader.getFileHeader();
			samHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
			
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
		} finally {
			reader.close();
		}
		
		log("Max read length is: " + readLength);
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
	
	private void sam2Fastq(String bam, String fastq, CompareToReference2 c2r, SAMFileWriter finalOutputSam) throws IOException {
		Sam2Fastq sam2Fastq = new Sam2Fastq();
		sam2Fastq.convert(bam, fastq, c2r, samHeader, finalOutputSam, this, regions);
//		if (isPairedEnd) {
//			sam2Fastq.convertPairedEnd(bam, fastq);
//		} else {
//			sam2Fastq.convert(bam, fastq);
//		}
	}
	
	private static int memCnt = 0;
	private void logOSMemory()  {
		
		System.out.println(memCnt++ +
				" mem total: " + Runtime.getRuntime().totalMemory() + 
				", max: " + Runtime.getRuntime().maxMemory() + 
				", free: " + Runtime.getRuntime().freeMemory());
	}
			
	private boolean cleanAndOutputContigs(String contigsSam, String cleanContigsFasta, boolean shouldRemoveSoftClips) throws IOException {
		
		Reference reference = new Reference(this.reference);
		
		boolean hasCleanContigs = false;
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(cleanContigsFasta, false));
		
		SAMFileReader contigReader = new SAMFileReader(new File(contigsSam));
		contigReader.setValidationStringency(ValidationStringency.SILENT);
		
		int contigCount = 0;
		
		for (SAMRecord contigRead : contigReader) {
			//TODO: Does this guarantee no alternate alignments?
//			if (contigRead.getMappingQuality() >= 1) {
			if (contigRead.getMappingQuality() >= this.minContigMapq) {
				
//				if (shouldRemoveSoftClips) {
					SAMRecordUtils.removeSoftClips(contigRead);
//				}
				
				String bases = contigRead.getReadString();
				
				//TODO: Why would cigar length be zero here?
				if ((bases.length() >= assemblerSettings.getMinContigLength()) &&
					(contigRead.getCigarLength() > 0)) {
					
					CigarElement first = contigRead.getCigar().getCigarElement(0);
					CigarElement last = contigRead.getCigar().getCigarElement(contigRead.getCigarLength()-1);
					
					if ((first.getOperator() == CigarOperator.M) &&
						(last.getOperator() == CigarOperator.M)) {

						// Pull in read length bases from reference to the beginning and end of the contig.
						String prefix = reference.getSequence(contigRead.getReferenceName(), 
								contigRead.getAlignmentStart()-readLength, readLength);
						String suffix = reference.getSequence(contigRead.getReferenceName(), contigRead.getAlignmentEnd()+1, readLength);
						
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
	
	// Assumes entirely soft clipped reads are filtered prior to here.
	// Checks only the first and last Cigar element
	// Does not adjust qualities
	
	private void alignToContigs(String tempDir, String inputSam, String alignedToContigSam,
			String contigFasta, CompareToReference2 c2r, SAMFileWriter finalOutputSam) throws IOException, InterruptedException {
		
		// Convert original bam to fastq
		String fastq = tempDir + "/" + "original_reads.fastq.gz";
		
		log("Preprocessing original reads for alignment: " + inputSam);
		System.out.println("c2r7: " + c2r);
		sam2Fastq(inputSam, fastq, c2r, finalOutputSam);
		log("Done preprocessing original reads for alignment: " + inputSam);
		
		//TODO: Only cut threads in half if in somatic mode.
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
	
	static class HitInfo {
		private SAMRecord record;
		private int position;
		private char strand;
		private int mismatches;
		
		public HitInfo(SAMRecord record, int position, char strand, int mismatches) {
			this.record = record;
			this.position = position;
			this.strand = strand;
			this.mismatches = mismatches;
		}

		public SAMRecord getRecord() {
			return record;
		}

		public int getPosition() {
			return position;
		}

		public char getStrand() {
			return strand;
		}
		
		public boolean isOnNegativeStrand() {
			return strand == '-';
		}
		
		public int getNumMismatches() {
			return mismatches;
		}
	}
	
	private int getIntAttribute(SAMRecord read, String attribute) {
		Integer num = read.getIntegerAttribute(attribute);
		
		if (num == null) {
			return 0;
		} else {
			return num;
		}
	}
	
	private int getOrigEditDistance(SAMRecord read) {
		
		Integer distance = null;
		
		if (read.getReadUnmappedFlag()) {
			distance = read.getReadLength();
		} else {
			distance = read.getIntegerAttribute("YX");
			
			if (distance == null) {
				distance = read.getReadLength();
			}
		}
				
		return distance;
	}
	
	public static int getEditDistance(SAMRecord read, CompareToReference2 c2r) {
		
		Integer distance = null;
		
		if (read.getReadUnmappedFlag()) {
			distance = read.getReadLength();
		} else if (c2r != null) {
			//distance = c2r.numMismatches(read) + getNumIndelBases(read);
			distance = c2r.numMismatches(read) + getNumIndelBases(read);
		} else {
			distance = read.getIntegerAttribute("NM");
			
			if (distance == null) {
				distance = read.getReadLength();
			}
		}
				
		return distance;
	}
	
	private RealignmentWriter getRealignmentWriter(SAMFileWriter outputReadsBam, boolean isTightAlignment, String tempDir) {
		RealignmentWriter writer;
		
		if (isTightAlignment && isPairedEnd) {
			writer = new BetaPairValidatingRealignmentWriter(this, outputReadsBam, tempDir);
		} else {
			writer = new SimpleRealignmentWriter(this, outputReadsBam, isTightAlignment);
		}
		
		return writer;
	}
	
	boolean isFiltered(SAMRecord read) {
		// Filter out single end reads when in paired end mode.
		return ((isPairedEnd) && (!read.getReadPairedFlag()));
	}
	
	private void spikeLog(String msg, SAMRecord read) {
		if (read.getReadName().contains("spike")) {
			System.out.println(msg);
		}
	}
	
	private int getClippingFactoredQuality(SAMRecord read) {
		
		// short circuit if cigar length is 1
		if (read.getCigarLength() == 1) {
			return read.getMappingQuality();
		}
		
		int mappedLength = read.getReadLength();

		CigarElement firstElement = read.getCigar().getCigarElement(0);
		CigarElement lastElement = read.getCigar().getCigarElement(read.getCigarLength()-1);
		
		if (firstElement.getOperator() == CigarOperator.S) {
			mappedLength -= firstElement.getLength();
		}
		
		if (lastElement.getOperator() == CigarOperator.S) {
			mappedLength -= lastElement.getLength();
		}
		
		int qual = read.getMappingQuality() * mappedLength / read.getReadLength();
		
		return qual;
	}
	
	private boolean isImprovedAlignment(SAMRecord read, SAMRecord orig, CompareToReference2 c2r) {
		
		boolean isImproved = false;
			
		Integer yr = orig.getIntegerAttribute("YR");
		if ((yr != null) && (yr == 1)) {
			// Original alignment was outside of target region list.
			// Calc updated edit distance to reference and compare to original
			
//			int origEditDistance = getOrigEditDistance(orig);
//			int updatedEditDistance = c2r.numMismatches(read) + getNumIndelBases(read);
			
			double origEditDistance = getOrigEditDistance(orig);
			double updatedEditDistance = c2r.numMismatches(read) + (1.5 * getNumIndels(read));
			
			isImproved = updatedEditDistance < origEditDistance;
		} else {
			isImproved = true;
		}
		
		return isImproved;
	}
	
	protected void adjustReads(String alignedToContigSam, SAMFileWriter outputSam, boolean isTightAlignment,
			CompareToReference2 c2r, String tempDir) throws IOException {
		
		log("Writing reads to: " + outputSam);
		
		RealignmentWriter writer = getRealignmentWriter(outputSam, isTightAlignment, tempDir);
		
		System.out.println("Opening sam file reader for: " + alignedToContigSam);
		System.out.flush();
		SAMFileReader contigReader = new SAMFileReader(new File(alignedToContigSam));
		System.out.println("Sam file reader opened for: " + alignedToContigSam);
		System.out.flush();
		contigReader.setValidationStringency(ValidationStringency.SILENT);
		
		SamStringReader samStringReader = new SamStringReader(samHeader);
		
		int ctr = 0;
		
		for (SAMRecord read : contigReader) {
			
			if ((ctr++ % 100000) == 0) {
				this.logOSMemory();
			}
			
			String origSamStr = read.getReadName();
			origSamStr = origSamStr.replace('|', '\t');
			SAMRecord orig;
			try {
				orig = samStringReader.getRead(origSamStr);
			} catch (RuntimeException e) {
				System.out.println("Error processing: [" + origSamStr + "]");
				System.out.println("Contig read: [" + read.getSAMString() + "]");
				e.printStackTrace();
				throw e;
			}
			orig.setHeader(samHeader);
			
			orig.setReadString(read.getReadString());
			orig.setBaseQualityString(read.getBaseQualityString());

			SAMRecord readToOutput = null;
			
//			if (orig.getReadUnmappedFlag()  && read.getReadName().contains("spike")) {
//				System.out.println("unmapped...");
//			}
			
			// Only adjust reads that align to contig with no indel and shorter edit distance than the original alignment
			String matchingString = read.getReadLength() + "M";
			if ((read.getCigarString().equals(matchingString)) &&
				(read.getReadUnmappedFlag() == false)  &&
				(!orig.getCigarString().contains("N")) &&  // Don't remap introns
				(getEditDistance(read, null) < getOrigEditDistance(orig)) &&
//				isImprovedAlignment(read, orig, c2r) &&
				(!isFiltered(orig))) {
				
				SAMRecord origRead = orig;
				String contigReadStr = read.getReferenceName();
				
				// Get alternate best hits
				List<HitInfo> bestHits = new ArrayList<HitInfo>();

				contigReadStr = contigReadStr.substring(contigReadStr.indexOf('~')+1);
				contigReadStr = contigReadStr.replace('~', '\t');
				SAMRecord contigRead = samStringReader.getRead(contigReadStr);
				
				int bestMismatches = getIntAttribute(read, "XM");
				
				// Filter this hit if it aligns past the end of the contig
				// Must use cigar length instead of read length, because the
				// the contig read bases are not loaded.
				if (read.getAlignmentEnd() <= contigRead.getCigar().getReadLength()) {
					HitInfo hit = new HitInfo(contigRead, read.getAlignmentStart(),
							read.getReadNegativeStrandFlag() ? '-' : '+', bestMismatches);
					
					bestHits.add(hit);
				}
				
				int numBestHits = getIntAttribute(read, "X0");
				int numSubOptimalHits = getIntAttribute(read, "X1");
				
				int totalHits = numBestHits + numSubOptimalHits;
				
				//TODO: If too many best hits, what to do?
				
//				spikeLog("total hits: " + totalHits, origRead);
				
				if ((totalHits > 1) && (totalHits < 1000)) {
					
//					if (totalHits < -1000) {
					// Look in XA tag.
					String alternateHitsStr = (String) read.getAttribute("XA");
					if (alternateHitsStr == null) {
						String msg = "best hits = " + numBestHits + ", but no XA entry for: " + read.getSAMString();
						System.out.println(msg);							
						this.missingXATag += 1;
					} else {
						
						String[] alternates = alternateHitsStr.split(";");
						
						for (int i=0; i<alternates.length-1; i++) {
							String[] altInfo = alternates[i].split(",");
							String altContigReadStr = altInfo[0];
							char strand = altInfo[1].charAt(0);
							int position = Integer.parseInt(altInfo[1].substring(1));
							String cigar = altInfo[2];
							int mismatches = Integer.parseInt(altInfo[3]);
							
							if ((cigar.equals(matchingString)) && (mismatches < bestMismatches)) {
								System.out.println("MISMATCH_ISSUE: " + read.getSAMString());
							}
							
							if ((cigar.equals(matchingString)) && (mismatches == bestMismatches)) {
								altContigReadStr = altContigReadStr.substring(altContigReadStr.indexOf('~')+1);
								altContigReadStr = altContigReadStr.replace('~', '\t');
								contigRead = samStringReader.getRead(altContigReadStr);
								
								// Filter this hit if it aligns past the end of the contig
								// Must use cigar length instead of read length, because the
								// the contig read bases are not loaded.
								if ((position + read.getReadLength()) <= contigRead.getCigar().getReadLength()) {
									HitInfo hit = new HitInfo(contigRead, position, strand, mismatches);
									bestHits.add(hit);
								}
							}
						}
					}
				}
				
				// chr_pos_cigar
				Map<String, SAMRecord> outputReadAlignmentInfo = new HashMap<String, SAMRecord>();
				
				for (HitInfo hitInfo : bestHits) {
					
					contigRead = hitInfo.getRecord();
					int position = hitInfo.getPosition() - 1;
					
					// Only consider this mapping if the assembled contig's quality is
					// greater than the original read's quality.
//					if (contigRead.getMappingQuality() > orig.getMappingQuality()) {
						
					//TODO: Normalize mapping quality
//					if (contigRead.getMappingQuality() > getClippingFactoredQuality(orig)) {
					
//					if (contigRead.getMappingQuality() > 30) {
					
					if (1 == 1) {

						List<ReadBlock> contigReadBlocks = ReadBlock.getReadBlocks(contigRead);
						
						ReadPosition readPosition = new ReadPosition(origRead, position, -1);
						SAMRecord updatedRead = updateReadAlignment(contigRead,
								contigReadBlocks, readPosition);
						
						if (updatedRead != null) {
							//TODO: Move into updateReadAlignment ?
//							if (updatedRead.getMappingQuality() == 0) {
//								updatedRead.setMappingQuality(1);
//							}
							
							if (updatedRead.getReadUnmappedFlag()) {
								updatedRead.setReadUnmappedFlag(false);
							}
							
							updatedRead.setReadNegativeStrandFlag(hitInfo.isOnNegativeStrand());

							// Set read bases to the aligned read (which will be expressed in 
							// forward strand context according to the primary alignment).
							// If this hit's strand is opposite the primary alignment, reverse the bases
//							if (hitInfo.isOnNegativeStrand() == read.getReadNegativeStrandFlag()) {
//								updatedRead.setReadString(read.getReadString());
//								updatedRead.setBaseQualityString(read.getBaseQualityString());
//							} else {
//								updatedRead.setReadString(reverseComplementor.reverseComplement(read.getReadString()));
//								updatedRead.setBaseQualityString(reverseComplementor.reverse(read.getBaseQualityString()));								
//							}
							
							// If the read's alignment info has been modified, record the original alignment.
							if (origRead.getReadUnmappedFlag() ||
								!origRead.getReferenceName().equals(updatedRead.getReferenceName()) ||
								origRead.getAlignmentStart() != updatedRead.getAlignmentStart() ||
								origRead.getReadNegativeStrandFlag() != updatedRead.getReadNegativeStrandFlag() ||
								!origRead.getCigarString().equals(updatedRead.getCigarString())) {
							
								if (isSoftClipEquivalent(origRead, updatedRead)) {
									// Restore Cigar and position
//									System.out.println("Re-setting [" + updatedRead.getSAMString() + "] --- [" + origRead.getSAMString() + "]");
									updatedRead.setAlignmentStart(origRead.getAlignmentStart());
									updatedRead.setCigar(origRead.getCigar());
									
								} else {
									String originalAlignment;
									if (origRead.getReadUnmappedFlag()) {
										originalAlignment = "N/A";
									} else {
										originalAlignment = origRead.getReferenceName() + ":" + origRead.getAlignmentStart() + ":" +
												(origRead.getReadNegativeStrandFlag() ? "-" : "+") + ":" + origRead.getCigarString();
									}
									
									// Read's original alignment position
									updatedRead.setAttribute("YO", originalAlignment);
								}
							}
							
							// Mismatches to the contig
							updatedRead.setAttribute("YM", hitInfo.getNumMismatches());
							
							// Contig's mapping quality
							updatedRead.setAttribute("YQ", hitInfo.getRecord().getMappingQuality());
							
							// Contig's length
							updatedRead.setAttribute("YL", hitInfo.getRecord().getCigar().getReadLength());
														
							// Check to see if this read has been output with the same alignment already.
//								String readAlignmentInfo;
							
							//TODO: Check strand!!!
							String readAlignmentInfo = updatedRead.getReferenceName() + "_" + updatedRead.getAlignmentStart() + "_" +
									updatedRead.getCigarString();
							/*
							if (isTightAlignment) {
								// Check ref, pos, strand, cigar
								readAlignmentInfo = updatedRead.getReferenceName() + "_" + updatedRead.getAlignmentStart() + "_" +
									(updatedRead.getReadNegativeStrandFlag() ? "-" : "+") + "_" + updatedRead.getCigarString();
							} else {
								// For loose alignment, allow read to align to either strand (pick one only)
								readAlignmentInfo = updatedRead.getReferenceName() + "_" + updatedRead.getAlignmentStart() + "_" +
										updatedRead.getCigarString();
							}
							*/
							
							if (!outputReadAlignmentInfo.containsKey(readAlignmentInfo)) {
								outputReadAlignmentInfo.put(readAlignmentInfo, updatedRead);
							}
						}
					}
				}
				
				if (outputReadAlignmentInfo.size() == 1) {
					//TODO: No need to iterate, just 1.  Dropping ambiguous reads
					readToOutput = outputReadAlignmentInfo.values().iterator().next();
					
					// Check to see if the original read location was in a non-target region.
					// If so, compare updated NM to ref versus original NM to ref
					if ((c2r != null) && (!isImprovedAlignment(readToOutput, orig, c2r))) {
						readToOutput = null;
					} else {
						int origBestHits = this.getIntAttribute(readToOutput, "X0");
						int origSuboptimalHits = this.getIntAttribute(readToOutput, "X1");
						
						// If the read mapped to multiple locations, set mapping quality to zero.
						if ((outputReadAlignmentInfo.size() > 1) || (totalHits > 1000)) {
							readToOutput.setMappingQuality(0);
						}
						
						// This must happen prior to updateMismatchAndEditDistance
						adjustForStrand(read.getReadNegativeStrandFlag(), readToOutput);
						
						if (readToOutput.getAttribute("YO") != null) {
							// HACK: Only add X0 for final alignment.  Assembler skips X0 > 1
							if (isTightAlignment) {
								readToOutput.setAttribute("X0", outputReadAlignmentInfo.size());
							} else {
								readToOutput.setAttribute("X0", null);
							}
							readToOutput.setAttribute("X1", origBestHits + origSuboptimalHits);
							
							// Clear various tags
							readToOutput.setAttribute("XO", null);
							readToOutput.setAttribute("XG", null);
							readToOutput.setAttribute("MD", null);
							readToOutput.setAttribute("XA", null);
							readToOutput.setAttribute("XT", null);
							
							if (c2r != null) {
								updateMismatchAndEditDistance(readToOutput, c2r, origRead);
							}
						}
					}
				}
			}
			
			adjustForStrand(read.getReadNegativeStrandFlag(), orig);
			writer.addAlignment(readToOutput, orig);
		}
//		System.out.println("Sleeping");
//		try {
//			Thread.sleep(200000);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}

		int realignedCount = writer.flush();
		contigReader.close();
		
		log("Done with: " + outputSam + ".  Number of reads realigned: " + realignedCount);
	}
	
	private boolean isSoftClipEquivalent(SAMRecord origRead, SAMRecord updatedRead) {
		
		boolean isEquivalent = false;
		
		if ((origRead.getCigarString().contains("S")) &&
			(origRead.getReferenceName().equals(updatedRead.getReferenceName())) &&
			(origRead.getReadNegativeStrandFlag() == updatedRead.getReadNegativeStrandFlag()) &&
			(origRead.getCigarLength() > 1)) {
			
			// Compare start positions
			int nonClippedOrigStart = origRead.getAlignmentStart();
			CigarElement firstElem = origRead.getCigar().getCigarElement(0); 
			if (firstElem.getOperator() == CigarOperator.S) {
				nonClippedOrigStart -= firstElem.getLength(); 
			}
			
			if (nonClippedOrigStart == updatedRead.getAlignmentStart()) {
				// Compare cigars
				List<CigarElement> elems = new ArrayList<CigarElement>(origRead.getCigar().getCigarElements());
				
				CigarElement first = elems.get(0);
				
				// If first element is soft clipped, lengthen the second element
				if (first.getOperator() == CigarOperator.S) {
					CigarElement second = elems.get(1);
					CigarElement newElem = new CigarElement(first.getLength() + second.getLength(), second.getOperator());
					elems.set(1,  newElem);
				}
				
				CigarElement last = elems.get(elems.size()-1);
				if (last.getOperator() == CigarOperator.S) {
					CigarElement nextToLast = elems.get(elems.size()-2);
					CigarElement newElem = new CigarElement(last.getLength() + nextToLast.getLength(), nextToLast.getOperator());
					elems.set(elems.size()-2, newElem);
				}
				
				List<CigarElement> newElems = new ArrayList<CigarElement>();

				for (CigarElement elem : elems) {
					if (elem.getOperator() != CigarOperator.S) {
						newElems.add(elem);
					}
				}
				
				Cigar convertedCigar = new Cigar(newElems);
				
				if (convertedCigar.equals(updatedRead.getCigar())) {
					isEquivalent = true;
				}
			}
		}
		
		return isEquivalent;
	}
	
	void adjustForStrand(boolean readAlreadyReversed, SAMRecord read) {
		if ( ((!readAlreadyReversed) && (read.getReadNegativeStrandFlag())) ||
			 ((readAlreadyReversed) && (!read.getReadNegativeStrandFlag())) ){
			read.setReadString(reverseComplementor.reverseComplement(read.getReadString()));
			read.setBaseQualityString(reverseComplementor.reverse(read.getBaseQualityString()));
		}
	}
	
	SAMRecord updateReadAlignment(SAMRecord contigRead,
			List<ReadBlock> contigReadBlocks, ReadPosition orig) {
		List<ReadBlock> blocks = new ArrayList<ReadBlock>();
		SAMRecord read = cloneRead(orig.getRead());

		read.setReferenceName(contigRead.getReferenceName());

		int contigPosition = orig.getPosition();
		int accumulatedLength = 0;

		// read block positions are one based
		// ReadPosition is zero based
		
		int totalInsertLength = 0;

		for (ReadBlock contigBlock : contigReadBlocks) {
			if ((contigBlock.getReadStart() + contigBlock.getReferenceLength()) >= orig
					.getPosition() + 1) {
				ReadBlock block = contigBlock.getSubBlock(accumulatedLength,
						contigPosition, read.getReadLength()
								- accumulatedLength);
				
				block.setReferenceStart(block.getReferenceStart() - totalInsertLength);
				
				// If this is an insert, we need to adjust the alignment start
				if ((block.getType() == CigarOperator.I) && (block.getLength() != 0)) {
					contigPosition = contigPosition - (contigBlock.getLength() - block.getLength());
					block.setReferenceStart(block.getReferenceStart() - (contigBlock.getLength() - block.getLength()));
//					block = contigBlock.getSubBlock(accumulatedLength,
//								contigPosition, read.getReadLength()
//								- accumulatedLength);
					
					totalInsertLength += block.getLength();
				}
				
				//TODO: Drop leading and trailing delete blocks

				// TODO: Investigate how this could happen
				if (block.getLength() != 0) {
					blocks.add(block);

					if (block.getType() != CigarOperator.D) {
						accumulatedLength += block.getLength();
					}

					if (accumulatedLength > read.getReadLength()) {
						throw new IllegalStateException("Accumulated Length: "
								+ accumulatedLength
								+ " is greater than read length: "
								+ read.getReadLength());
					}

					if (accumulatedLength == read.getReadLength()) {
						break;
					}
				}
			}
		}

		if (blocks.size() > 0) {
			
			// Don't allow reads to align past end of contig.
//			int cigarLength = ReadBlock.getTotalLength(blocks);
//			if (cigarLength != read.getReadLength()) {
//				read = null;
//			}
			
			// If we've aligned past the end of the contig resulting in a short Cigar
			// length, append additonal M to the Cigar			
			ReadBlock.fillToLength(blocks, read.getReadLength());
			int newAlignmentStart = blocks.get(0).getReferenceStart();
			String newCigar = ReadBlock.toCigarString(blocks);

			read.setCigarString(newCigar);
			read.setAlignmentStart(newAlignmentStart);
		} else {
			// TODO: Investigate how this could happen.
			read = null;
		}

		return read;
	}

	private SAMRecord cloneRead(SAMRecord read) {
		try {
			return (SAMRecord) read.clone();
		} catch (CloneNotSupportedException e) {
			// Infamous "this should never happen" comment here.
			e.printStackTrace();
			throw new RuntimeException(e);
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
					(currentRegion.getEnd() + readLength >= region.getStart())) {
					
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
			
			if (end < region.getEnd()) {
				end += regionOverlap;
			}
			
			// If we're at or near the end of the region, stop at region end.
			if (end > (region.getEnd() - minRegionRemainder)) {
				end = region.getEnd();
			}
			
			pos = marker;
			
			regions.add(new Feature(region.getSeqname(), start, end));
		}
		
		return regions;
	}
		
	private Assembler newAssembler() {
		//Assembler assem = new JavaAssembler();
		NativeAssembler assem = new NativeAssembler();
		
//		if (assem instanceof JavaAssembler) {
//			JavaAssembler ja = (JavaAssembler) assem;
//			ja.setKmerSize(assemblerSettings.getKmerSize());
//			ja.setMinNodeFrequncy(assemblerSettings.getMinNodeFrequncy());
//			ja.setMinContigLength(assemblerSettings.getMinContigLength());
//			ja.setMinContigRatio(assemblerSettings.getMinContigRatio());
//		}

		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(assemblerSettings
				.getMaxPotentialContigs());

		assem.setMaxPathsFromRoot(100000);
		assem.setReadLength(readLength);
		assem.setKmer(assemblerSettings.getKmerSize());
		assem.setMinKmerFrequency(assemblerSettings.getMinNodeFrequncy());

		return assem;
	}
	
	private Assembler newUnalignedAssembler(int mnfMultiplier) {
		//Assembler assem = new JavaAssembler();
		NativeAssembler assem = new NativeAssembler();
		
//		if (assem instanceof JavaAssembler) {
//			JavaAssembler ja = (JavaAssembler) assem;
//			ja.setKmerSize(assemblerSettings.getKmerSize());
//			ja.setMinNodeFrequncy(assemblerSettings.getMinUnalignedNodeFrequency() * mnfMultiplier);
//			ja.setMinContigLength(assemblerSettings.getMinContigLength());
//			ja.setMinContigRatio(-1.0);
//		}

		assem.setMaxContigs(MAX_POTENTIAL_UNALIGNED_CONTIGS);
		assem.setTruncateOutputOnRepeat(false);
		assem.setMaxPathsFromRoot(5000);
		assem.setReadLength(readLength);
		// Could be smaller for higher sensitivity here?
		assem.setKmer(assemblerSettings.getKmerSize());
		assem.setMinKmerFrequency(assemblerSettings.getMinUnalignedNodeFrequency());

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
	
	public SAMFileHeader getHeader() {
		return this.samHeader;
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

	public static void run(String[] args) throws Exception {
		
		System.out.println("Starting 0.51 ...");
		
		ReAlignerOptions options = new ReAlignerOptions();
		options.parseOptions(args);

		if (options.isValid()) {

			AssemblerSettings assemblerSettings = new AssemblerSettings();

			assemblerSettings.setKmerSize(options.getKmerSize());
			assemblerSettings.setMinContigLength(options.getMinContigLength());
			assemblerSettings.setMinNodeFrequncy(options.getMinNodeFrequency());
			assemblerSettings.setMaxPotentialContigs(options
					.getMaxPotentialContigs());
			assemblerSettings.setMinContigRatio(options.getMinContigRatio());
			assemblerSettings.setMinUnalignedNodeFrequency(options.getMinUnalignedNodeFrequency());

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

			realigner.reAlign(options.getInputFile(), options.getInputFile2(), options.getInputFile3(),
					options.getOutputFile(), options.getOutputFile2(), options.getOutputFile3());

			long e = System.currentTimeMillis();

			System.out.println("Elapsed seconds: " + (e - s) / 1000);
		}
	}
	
	public static void main(String[] args) throws Exception {
		ReAligner realigner = new ReAligner();
//		String originalReadsSam = args[0];
//		String alignedToContigSam = args[1];
//		String unalignedSam = args[2];
//		String outputFilename = args[3];
		
		String originalReadsSam = "/home/lmose/dev/ayc/sim/s43/orig_atc.bam";
		String alignedToContigSam = "/home/lmose/dev/ayc/sim/s43/atc.bam";
//		String unalignedSam = args[2];
		String outputFilename = "/home/lmose/dev/ayc/sim/s43/atc_out.bam";
		
		SAMFileReader reader = new SAMFileReader(new File(originalReadsSam));
	
		realigner.samHeader = reader.getFileHeader();
		
		reader.close();
		
		realigner.readLength = 76;
		realigner.isPairedEnd = true;
		realigner.minInsertLength = 70;
		realigner.maxInsertLength = 550;
		realigner.regionsGtf = "/home/lmose/dev/ayc/p3/wxs.gtf";
		realigner.loadRegions();
		
//		realigner.getSamHeader(originalReadsSam);
		
		SAMFileWriter outputReadsBam = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				realigner.samHeader, true, new File(outputFilename));
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("/home/lmose/reference/chr19/19.fa");
		
		SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				realigner.samHeader, false, new File("/home/lmose/dev/ayc/p3/output.bam"));

		realigner.adjustReads("/home/lmose/dev/ayc/p3/a2c2.bam", writer, true, c2r, "/home/lmose/dev/ayc/p3");
		
		outputReadsBam.close();
	}


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

