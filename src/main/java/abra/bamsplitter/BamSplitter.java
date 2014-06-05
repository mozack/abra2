package abra.bamsplitter;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import abra.ThreadManager;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

public class BamSplitter {

	public void split(String filename, int numThreads, String outputDirectory) throws IOException, InterruptedException {
		long s = System.currentTimeMillis();
		
		File dir = new File(outputDirectory);
		if (!dir.exists()) {
			dir.mkdir();
		}
		
		SAMFileReader rdr = new SAMFileReader(new File(filename));
		
		ThreadManager threads = new ThreadManager(numThreads);
		
		Map<String, SAMFileWriter> outputWriterMap = new HashMap<String, SAMFileWriter>();
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		writerFactory.setUseAsyncIo(false);
		
		// Farm each chromosome out to its own thread.
		for (SAMSequenceRecord chr : rdr.getFileHeader().getSequenceDictionary().getSequences()) {		
			SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
					rdr.getFileHeader(), false, new File(outputDirectory + "/" + chr.getSequenceName() + ".bam"));
			
			outputWriterMap.put(chr.getSequenceName(), writer);
			
			BamSplitterThread thread = new BamSplitterThread(threads, filename, chr.getSequenceName(), writer);
			threads.spawnThread(thread);
		}
		threads.waitForAllThreadsToComplete();
		
		// Now go back and retrieve the unmapped reads.
		System.err.println("Processing unmapped reads");
		Iterator<SAMRecord> iter = rdr.queryUnmapped();
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
		
			// If this read is not assigned a position, but the mate is, include in the output BAM associated with mate's chromosome.
			if (read.getReferenceIndex() == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && read.getMateReferenceIndex() != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
				SAMFileWriter writer = outputWriterMap.get(read.getMateReferenceName());
				writer.addAlignment(read);
			}
		}
		
		for (SAMFileWriter writer : outputWriterMap.values()) {
			writer.close();
		}
		
		rdr.close();
		
		long e = System.currentTimeMillis();
		
		System.err.println("BAMSplitter done.  Elapsed minutes: " + (double) (e-s)/1000.0/60.0);
	}
	
	public static void main(String[] args) throws Exception {
//		int numThreads = Integer.parseInt(args[0]);
//		String inputFile = args[1];
//		String outputDir = args[2];

		int numThreads = 2;
		String inputFile = "/home/lmose/dev/abra/splitter/tumor.sort.bam";
		String outputDir = "/home/lmose/dev/abra/splitter/split";
		
		new BamSplitter().split(inputFile, numThreads, outputDir);
	}
}
