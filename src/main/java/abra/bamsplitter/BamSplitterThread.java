package abra.bamsplitter;

import java.io.File;
import java.util.Iterator;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import abra.AbraRunnable;
import abra.ThreadManager;

public class BamSplitterThread extends AbraRunnable {
	
	private String inputFile;
	private String chromosome;
	private SAMFileWriter outputWriter;

	public BamSplitterThread(ThreadManager threadManager, String inputFile, String chromosome, SAMFileWriter outputWriter) {
		super(threadManager);
		
		this.inputFile = inputFile;
		this.outputWriter = outputWriter;
		this.chromosome = chromosome;
	}
	
	@Override
	public void go() throws Exception {
		System.err.println("Starting chromosome: " + chromosome);
		SAMFileReader rdr = new SAMFileReader(new File(inputFile));
		rdr.setValidationStringency(ValidationStringency.SILENT);
		
		Iterator<SAMRecord> iter = rdr.queryContained(chromosome, 0, 0);
			
		while (iter.hasNext()) {
			SAMRecord read = (SAMRecord) iter.next();
			
			outputWriter.addAlignment(read);
		}
		
		rdr.close();
		System.err.println("Chromosome: " + chromosome + " done.");
	}
}
