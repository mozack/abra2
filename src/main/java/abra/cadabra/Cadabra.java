package abra.cadabra;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import abra.CompareToReference2;
import abra.Feature;
import abra.Logger;
import abra.SAMRecordUtils;
import abra.ThreadManager;

import abra.cadabra.CadabraProcessor.SampleCall;
import abra.cadabra.CadabraProcessor.SomaticCall;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;

public class Cadabra {

	private CompareToReference2 c2r;
	
	private Map<String, List<SampleCall>> chromosomeCalls = new HashMap<String, List<SampleCall>>();
	private Map<String, List<SomaticCall>> chromosomeSomaticCalls = new HashMap<String, List<SomaticCall>>();
	
	public void call(CadabraOptions options) throws IOException, InterruptedException {
		c2r = new CompareToReference2();
		c2r.init(options.getReference());
		
		outputHeader(options);
		
		ThreadManager threadManager = new ThreadManager(options.getNumThreads());
		
		for (String chromosome : c2r.getChromosomes()) {
			Feature region = new Feature(chromosome, 1, c2r.getChromosomeLength(chromosome));
			CadabraRunnable thread = new CadabraRunnable(threadManager, this, options, c2r, region);			
			threadManager.spawnThread(thread);
		}
		
		threadManager.waitForAllThreadsToComplete();
		
		// Output calls.
		if (options.getNormal() == null) {
			// Simple calling
			for (String chromosome : c2r.getChromosomes()) {
				for (SampleCall call : chromosomeCalls.get(chromosome)) {
					System.out.println(call);
				}
			}
		} else {
			// Somatic calling
			for (String chromosome : c2r.getChromosomes()) {
				for (SomaticCall call : chromosomeSomaticCalls.get(chromosome)) {
					System.out.println(call);
				}
			}
		}
		
		Logger.info("Cadabra done.");
	}
	
	void addCalls(String chromosome, List<SampleCall> calls) {
		Logger.info("Choromosome: %s done.", chromosome);
		synchronized(chromosomeCalls) {
			chromosomeCalls.put(chromosome, calls);
		}
	}
	
	void addSomaticCalls(String chromosome, List<SomaticCall> calls) {
		Logger.info("Choromosome: %s done.", chromosome);
		synchronized(chromosomeSomaticCalls) {
			chromosomeSomaticCalls.put(chromosome, calls);
		}
	}
	
	private void outputHeader(CadabraOptions options) throws IOException {
		
		SAMFileHeader header;
		String vcfColumns;
		
		if (options.getNormal() == null) {
			// Single sample
			SamReader reader = SAMRecordUtils.getSamReader(options.getTumor());
			header = reader.getFileHeader();
			reader.close();
			vcfColumns = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE";
		} else {
			// Somatic
			SamReader reader = SAMRecordUtils.getSamReader(options.getNormal());
			SAMFileHeader normalHeader = reader.getFileHeader();
			reader.close();
			header = normalHeader;
			
			reader = SAMRecordUtils.getSamReader(options.getTumor());
			SAMFileHeader tumorHeader = reader.getFileHeader();
			reader.close();
			
			//TODO: Double check against specified reference?
			if (!normalHeader.getSequenceDictionary().equals(tumorHeader.getSequenceDictionary())) {
				Logger.error("Reference Sequences for tumor and normal do not match.  Check the VCF headers.");
			}
			
			vcfColumns = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR";
		}
		
		System.out.println("##fileformat=VCFv4.2");
		System.out.println("##reference=file://" + options.getReference());
		
		for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
			System.out.println(String.format("##contig=<ID=%s,length=%d>", seq.getSequenceName(), seq.getSequenceLength()));
		}
		
		System.out.println("##INFO=<ID=RP,Number=1,Type=Integer,Description=\"Number of times smallest repeating alternate sequence appears in the reference\">");
		System.out.println("##INFO=<ID=RU,Number=1,Type=String,Description=\"Smallest repeat unit within alternate sequence.  Appears RP times in reference\">");
		System.out.println("##INFO=<ID=HRUN,Number=2,Type=Integer,Description=\"Length,position of homopolymer run found in CTX\">");
		System.out.println("##INFO=<ID=CTX,Number=1,Type=String,Description=\"Reference context sequence\">");
		System.out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		System.out.println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth (fragment)\">");
		System.out.println("##FORMAT=<ID=DP2,Number=1,Type=Integer,Description=\"Depth 2 (read)\">");
		System.out.println("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele Depth (fragment)\">");
		System.out.println("##FORMAT=<ID=AD2,Number=2,Type=Integer,Description=\"Allele Depth (read)\">");
		System.out.println("##FORMAT=<ID=ROR,Number=4,Type=Integer,Description=\"Read Orientation (ref_fwd, ref_rev, alt_fwd, alt_rev)\">");
		System.out.println("##FORMAT=<ID=LMQ,Number=1,Type=Integer,Description=\"Number of reads filtered due to low mapping quality\">");
		System.out.println("##FORMAT=<ID=ISPAN,Number=1,Type=Integer,Description=\"Max variant read pos minus min variant read pos\">");
		System.out.println("##FORMAT=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele frequency\">");
		System.out.println("##FORMAT=<ID=MER,Number=1,Type=Integer,Description=\"Number of ref reads with num mismatches greater than read length * .05\">");
		System.out.println("##FORMAT=<ID=FROR,Number=1,Type=Float,Description=\"Phred scaled Fisher's Exact Test for read orientation\">");
		System.out.println(vcfColumns);
	}
		
	public static void main(String[] args) throws Exception {
//		String normal = "/home/lmose/dev/abra/cadabra/normal_test2.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/tumor_test2.bam";
		
//		String reference = "/home/lmose/reference/chr1/1.fa";
//		String normal = "/home/lmose/dev/abra/cadabra/normal1.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/tumor1.bam";

		
//		String normal = "/home/lmose/dev/abra/cadabra/normal.abra4.sort.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/tumor.abra4.sort.bam";

//		String reference = "/home/lmose/reference/chr1/chr1.fa";
//		String normal = "/home/lmose/dev/abra/cadabra/t2/ntest.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/t2/ttest.bam";

		
//		String reference = "/home/lmose/reference/chr1/chr1.fa";
//		String normal = "/home/lmose/dev/abra/cadabra/ins/ntest.bam";
//		String tumor = "/home/lmose/dev/abra/cadabra/ins/ttest.bam";
		
		CadabraOptions options = new CadabraOptions();
		options.parseOptions(args);
		
		if (options.isValid()) {
			new Cadabra().call(options);
		}
	}
}
