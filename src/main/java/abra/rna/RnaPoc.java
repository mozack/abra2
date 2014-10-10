package abra.rna;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import abra.NativeAssembler;
import abra.NativeLibraryLoader;
import abra.ReverseComplementor;
import abra.SAMRecordUtils;
import abra.ThreadManager;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class RnaPoc {
	
	//public static int MAX_READ_GAP = 1000000;
	public static int MAX_READ_GAP = 50000;
	
	private BufferedWriter contigWriter;
	private BufferedWriter badRegionBed;
	
//	private BufferedWriter reads1;
//	private BufferedWriter reads2;
	
	private ThreadManager threadManager;
	
	private void init(String tempDir) {
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

	public void run(String input, String output, String temp, int numThreads, String badRegionBedFile) throws IOException, InterruptedException {
		
		init(temp);
		
		this.threadManager = new ThreadManager(numThreads);
		
		contigWriter = new BufferedWriter(new FileWriter(output, false));
		badRegionBed = new BufferedWriter(new FileWriter(badRegionBedFile, false));
		
//		reads1 = new BufferedWriter(new FileWriter(readsFile + "1.fa", false));
//		reads2 = new BufferedWriter(new FileWriter(readsFile + "2.fa", false));
		
		List<SAMRecord> currReads = new ArrayList<SAMRecord>();
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		int prevMaxEnd = -1;
		
		SAMRecord lastRead = null;

		for (SAMRecord read : reader) {
			if (read.getMappingQuality() > 0) {
				if (lastRead == null || !lastRead.getReferenceName().equals(read.getReferenceName()) || (read.getAlignmentStart()-prevMaxEnd) < MAX_READ_GAP) {
					currReads.add(read);
				} else {
//					processReads(currReads);
					spawnProcessingThread(currReads);
					currReads = new ArrayList<SAMRecord>();
					currReads.add(read);
				}
				
				if (read.getAlignmentEnd() > prevMaxEnd || !lastRead.getReferenceName().equals(read.getReferenceName())) {
					prevMaxEnd = read.getAlignmentEnd();
				}
				
				lastRead = read;
			}
		}
		
		threadManager.waitForAllThreadsToComplete();
		
		reader.close();
		contigWriter.close();
		badRegionBed.close();
	}
	
	private void spawnProcessingThread(List<SAMRecord> reads) {
		RnaRegionHandler handler = new RnaRegionHandler(threadManager, this, reads);
		threadManager.spawnThread(handler);
	}
	
	private String getRegionEntry(List<SAMRecord> reads) {
		String chr = reads.get(0).getReferenceName();
		int start = reads.get(0).getAlignmentStart();
		
		int stop = start + 1;
		for (SAMRecord read : reads) {
			if (read.getAlignmentEnd() > stop) {
				stop = read.getAlignmentEnd();
			}
		}
		
		return chr + "\t" + start + "\t" + stop;
	}
	
	private void filterNbases(List<SAMRecord> reads) {
		Iterator<SAMRecord> iter = reads.iterator();
		
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			if (read.getReadString().contains("N")) {
				iter.remove();
			}
		}
		
	}
	
	void processReads(List<SAMRecord> reads) throws IOException {

		filterNbases(reads);
		
		if (reads.size() > 1) {
			NativeAssembler assem = newAssembler();
			
			String contigs = assem.simpleAssemble(reads);
			
			if (contigs.equals("<ERROR>") || contigs.equals("<REPEAT>")) {
				
				badRegionBed.write(getRegionEntry(reads) + "\n");
				
				/*
				// Pair and output original reads
				Collections.sort(reads, new ReadNameComparator());
				
				String name = "";
				SAMRecord first = null;
				SAMRecord second = null;
				
				for (SAMRecord read : reads) {
					if (SAMRecordUtils.isPrimary(read)) {
						if (read.getFirstOfPairFlag()) {
							first = read;
							if (second != null) {
								if (first.getReadName().equals(second.getReadName())) {
									outputRead(reads1, first);
									outputRead(reads2, second);
									first = null;
									second = null;
								} else {
									second = null;
								}
							}
						} else {
							second = read;
							if (first != null) {
								if (first.getReadName().equals(second.getReadName())) {
									outputRead(reads1, first);
									outputRead(reads2, second);
									first = null;
									second = null;
								} else {
									first = null;
								}
							}
						}
					}
				}
				*/
				
			} else if (!contigs.isEmpty()) {
				appendContigs(contigs);
			}
		}
	}
	
	private void outputRead(BufferedWriter writer, SAMRecord read) throws IOException {
		writer.write(">" + read.getReadName() + "\n");
		String bases = read.getReadString();
		if (read.getReadNegativeStrandFlag()) {
			bases = new ReverseComplementor().reverseComplement(bases);
		}
		writer.write(bases + "\n");
	}
	
	private synchronized void appendContigs(String contigs) throws IOException {
		contigWriter.write(contigs);
	}
	
	private NativeAssembler newAssembler() {
		NativeAssembler assem = new NativeAssembler();

		assem.setTruncateOutputOnRepeat(true);
		assem.setMaxContigs(10000);

		assem.setMaxPathsFromRoot(5000000);
		assem.setReadLength(75);
		assem.setKmer(new int[] { 17, 27, 37, 47 });
		assem.setMinKmerFrequency(2);
		assem.setMinBaseQuality(40);
		
		// The following params not used
		assem.setMinReadCandidateFraction(0);
		assem.setMaxAverageDepth(0);
		assem.setShouldSearchForSv(false);
		assem.setAverageDepthCeiling(0);

		return assem;		
	}
	
	static class ReadNameComparator implements Comparator<SAMRecord> {

		@Override
		public int compare(SAMRecord o1, SAMRecord o2) {
			return o1.getReadName().compareTo(o2.getReadName());
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		RnaPoc poc = new RnaPoc();
		String input = args[0];
		String output = args[1];
		String temp = args[2];
		int numThreads = Integer.parseInt(args[3]);
		String readsFile = args[4];
		
		poc.run(input, output, temp, numThreads, readsFile);
	}
}
