package abra;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class NativeAssembler implements Assembler {
	
	private boolean truncateOnRepeat;
	private int maxContigs;
	private int maxPathsFromRoot;
	private int readLength;
	private int kmer;
	private Set<String> readIds;

	private native int assemble(String input, String output, String prefix, int truncateOnRepeat, int maxContigs, int maxPathsFromRoot, int readLength, int kmerSize);
	
	static{
        System.loadLibrary("Abra");
	}
	
	private String getIdentifier(SAMRecord read) {
		String id = read.getReadName();
		
		if (read.getReadPairedFlag() && read.getSecondOfPairFlag()) {
			id += "_2";
		}
		
		return id;
	}
	
//	public boolean assembleContigs(String input, String output, String prefix, boolean checkForDupes) {
	
	public boolean assembleContigs(List<String> inputFiles, String output, String tempDir, Feature region, String prefix, boolean checkForDupes, ReAligner realigner) {
		
		long start = System.currentTimeMillis();
		
		int count = 0;
		
		readIds = new HashSet<String>();
		
		try {
			
			String readFile;
			if (region != null) {
				readFile = tempDir + "/" + region.getDescriptor() + ".reads";
			} else {
				readFile = tempDir + "/" + "unaligned.reads";
			}
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(readFile, false));

			for (String input : inputFiles) {
				SAMFileReader reader = new SAMFileReader(new File(input));
				reader.setValidationStringency(ValidationStringency.SILENT);
	
				Iterator<SAMRecord> iter;
				if (region != null) {
					iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
				} else {
					iter = reader.iterator();
				}
				
				while (iter.hasNext()) {
					
					SAMRecord read = iter.next();
					
					if (read.getReadLength() > readLength) {
						reader.close();
						throw new IllegalArgumentException(
								"Read length exceeds expected value of: " + readLength + " for read [" +
								read.getSAMString() + "]");
					}
					
					// Don't allow same read to be counted twice.
					if ( (!realigner.isFiltered(read)) && (!read.getDuplicateReadFlag()) && (!read.getReadFailsVendorQualityCheckFlag()) && ((!checkForDupes) || (!readIds.contains(getIdentifier(read))))) {
	//					boolean hasAmbiguousBases = read.getReadString().contains("N");
						Integer numBestHits = (Integer) read.getIntegerAttribute("X0");
						boolean hasAmbiguousInitialAlignment = numBestHits != null && numBestHits > 1;
						//TODO: Stampy ambiguous read (mapq < 4)
						
	//					if (!hasAmbiguousBases && !hasAmbiguousInitialAlignment && !hasLowQualityBase(read)) {
						if (!hasAmbiguousInitialAlignment) {
							if (!checkForDupes) {
								readIds.add(getIdentifier(read));
							}
							
							if (read.getReadLength() == readLength) {
								writer.write(read.getReadString() + "\n");
								writer.write(read.getBaseQualityString() + "\n");
							} else {
								StringBuffer basePadding = new StringBuffer();
								StringBuffer qualPadding = new StringBuffer();
								
								for (int i=0; i<readLength-read.getReadLength(); i++) {
									basePadding.append('N');
									qualPadding.append('!');
								}
								
								writer.write(read.getReadString() + basePadding.toString() + "\n");
								writer.write(read.getBaseQualityString() + qualPadding.toString() + "\n");							
							}
						}
					}
				}
				
				reader.close();
			}
			
			readIds = null;
			writer.close();
			
			count = assemble(
					readFile,
					output, 
					prefix, 
					truncateOnRepeat ? 1 : 0,
					maxContigs,
					maxPathsFromRoot,
					readLength,
					kmer);
			
			File inputReadFile = new File(readFile);
			inputReadFile.delete();

		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
		
		long end = System.currentTimeMillis();
		
		System.out.println("Elapsed msecs in NativeAssembler: " + (end-start));
		
		return count != 0;
	}
	
	private boolean hasLowQualityBase(SAMRecord read) {
		//TODO: Don't hardcode phred33
		for (int i=0; i<read.getBaseQualityString().length(); i++) {
			if ((read.getBaseQualityString().charAt(i) - '!') < 20) {
				return true;
			}
		}
		
		return false;
	}

	public boolean isTruncateOnRepeat() {
		return truncateOnRepeat;
	}

	public void setTruncateOutputOnRepeat(boolean truncateOnRepeat) {
		this.truncateOnRepeat = truncateOnRepeat;
	}

	public int getMaxContigs() {
		return maxContigs;
	}

	public void setMaxContigs(int maxContigs) {
		this.maxContigs = maxContigs;
	}

	public int getMaxPathsFromRoot() {
		return maxPathsFromRoot;
	}

	public void setMaxPathsFromRoot(int maxPathsFromRoot) {
		this.maxPathsFromRoot = maxPathsFromRoot;
	}
	
	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}
	
	public void setKmer(int kmer) {
		this.kmer = kmer;
	}
	
//	public static void run(String input, String output) {
//		NativeAssembler assem = new NativeAssembler();
//		assem.setTruncateOutputOnRepeat(false);
//		assem.setMaxContigs(2000000);
//		assem.setMaxPathsFromRoot(5000);
//		
//		assem.assembleContigs(input, output, "contig", true);
//		
//	}
	
	public static void main(String[] args) {
		NativeAssembler assem = new NativeAssembler();
		assem.setTruncateOutputOnRepeat(false);
		assem.setMaxContigs(2000000);
		assem.setMaxPathsFromRoot(5000);
		/*
		assem.assembleContigs(args[0], args[1], "contig");
		*/
		
//		for (int i=0; i<10; i++) {
//			run(args[0], args[1] + "_" + i);
//		}
		
//		run(args[0], args[1]);
		
//		assem.assembleContigs("/home/lmose/code/abra/src/main/c/1810_reads.txt",
//				"/home/lmose/code/abra/src/main/c/1810.fa", "bar");
	}
}
