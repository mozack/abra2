package abra.cadabra;

import abra.Logger;
import abra.Options;
import joptsimple.OptionParser;

public class CadabraOptions extends Options {
	
	private static final String NUM_THREADS = "threads";
	private static final String TUMOR = "tumor";
	private static final String NORMAL = "normal";
	private static final String SAMPLE = "sample";
	private static final String REFERENCE = "ref";
	private static final String STRP_THRESHOLD = "strp";
	private static final String HRUN_THRESHOLD = "hrun";
	private static final String ISPAN_FILTER = "ispan";
	// Filter variants below this qual threshold
	private static final String QUAL_FILTER = "qual";
	private static final String FS_FILTER = "fs";
	private static final String LOW_MQ_FILTER = "lmq";
	// Do not output variants below this qual threshold
	private static final String MIN_QUAL = "mq";
	private static final String MIN_MAPQ = "mapq";
	private static final String MIN_VAF = "mf";
	private static final String PCR_PENALTY = "pen";
	

	private OptionParser parser;
	private boolean isValid;
	
	private int numThreads;
	private String tumor;
	private String normal;
	private String reference;
	private int strpThreshold;
	private int hrunThreshold;
	private int pcrPenalty;
	private int ispanFilter;
	private float qualFilter;
	private int fsFilter;
	private float lowMQFilter;
	private float minQual;
	private float minMapq;
	private float minVaf;
	
	@Override
	protected OptionParser getOptionParser() {
    	if (parser == null) {
            parser = new OptionParser();
            parser.accepts(NUM_THREADS, "Processing threads").withRequiredArg().ofType(Integer.class).defaultsTo(4);
            parser.accepts(TUMOR, "Tumor BAM file").withRequiredArg().ofType(String.class);
            parser.accepts(NORMAL, "Normal BAM file").withRequiredArg().ofType(String.class);
            parser.accepts(SAMPLE, "BAM file to be used for single sample calling").withRequiredArg().ofType(String.class);
            parser.accepts(REFERENCE, "Reference fasta").withRequiredArg().ofType(String.class);
            parser.accepts(STRP_THRESHOLD, "Filter variants with short tandem repeat period at or above this threshold(-1 to disable)").withRequiredArg().ofType(Integer.class).defaultsTo(5); 
            parser.accepts(HRUN_THRESHOLD, "Filter short indels with a nearby homopolymer run of this length or greater - only neighboring 20 bases searched (-1 to disable)").withRequiredArg().ofType(Integer.class).defaultsTo(6);
            parser.accepts(PCR_PENALTY, "Penalize quality score for variants reaching strp or hrun thresholds by specified amount.").withRequiredArg().ofType(Integer.class).defaultsTo(30);
            parser.accepts(ISPAN_FILTER, "Filter variants with max index span less than specified value").withRequiredArg().ofType(Integer.class).defaultsTo(19);
            parser.accepts(QUAL_FILTER, "Filter variants with quality score less than specified value").withRequiredArg().ofType(Float.class).defaultsTo(5f);
            parser.accepts(FS_FILTER, "Filter variants with FS score greater than specified value").withRequiredArg().ofType(Integer.class).defaultsTo(70);
            parser.accepts(LOW_MQ_FILTER, "Filter variants with fraction of low quality reads greater than specified value").withRequiredArg().ofType(Float.class).defaultsTo(.5f);
            parser.accepts(MIN_QUAL, "Variants with quality below specified threshold are not output").withRequiredArg().ofType(Float.class).defaultsTo(5.0f);
            parser.accepts(MIN_MAPQ, "Reads with mapping quality below specified value are excluded from processing (except unmapped reads)").withRequiredArg().ofType(Integer.class).defaultsTo(20);
            parser.accepts(MIN_VAF, "Do not output variants with frequency below specified value").withRequiredArg().ofType(Float.class).defaultsTo(0.0f);            
    	}
    	
    	return parser;
	}
	
	public void init() {
		this.numThreads = (Integer) getOptions().valueOf(NUM_THREADS);
		this.tumor = (String) getOptions().valueOf(TUMOR);
		if (tumor == null) {
			this.tumor = (String) getOptions().valueOf(SAMPLE);
		}
		this.normal = (String) getOptions().valueOf(NORMAL);
		this.reference = (String) getOptions().valueOf(REFERENCE);
		this.strpThreshold = (Integer) getOptions().valueOf(STRP_THRESHOLD);
		this.hrunThreshold = (Integer) getOptions().valueOf(HRUN_THRESHOLD);
		this.ispanFilter = (Integer) getOptions().valueOf(ISPAN_FILTER);
		this.qualFilter = (Float) getOptions().valueOf(QUAL_FILTER);
		this.fsFilter = (Integer) getOptions().valueOf(FS_FILTER);
		this.lowMQFilter = (Float) getOptions().valueOf(LOW_MQ_FILTER);
		this.minQual = (Float) getOptions().valueOf(MIN_QUAL);
		this.minMapq = (Integer) getOptions().valueOf(MIN_MAPQ);
		this.minVaf = (Float) getOptions().valueOf(MIN_VAF);
		this.pcrPenalty = (Integer) getOptions().valueOf(PCR_PENALTY);
	}

	@Override
	protected void validate() {
		
		isValid = true;
		
		if (tumor == null) {
			Logger.error("Please specify a " + TUMOR + " or a " + SAMPLE);
			isValid = false;
		}
		
		if (reference == null) {
			Logger.error("Please specify a reference");
			isValid = false;
		}
	}
	
	public boolean isValid() {
		return isValid;
	}

	public int getNumThreads() {
		return numThreads;
	}

	public String getTumor() {
		return tumor;
	}

	public String getNormal() {
		return normal;
	}

	public String getReference() {
		return reference;
	}

	public int getStrpThreshold() {
		return strpThreshold;
	}

	public int getHrunThreshold() {
		return hrunThreshold;
	}
	
	public int getPcrPenalty() {
		return pcrPenalty;
	}

	public int getIspanFilter() {
		return ispanFilter;
	}

	public float getQualFilter() {
		return qualFilter;
	}

	public int getFsFilter() {
		return fsFilter;
	}

	public float getLowMQFilter() {
		return lowMQFilter;
	}

	public float getMinQual() {
		return minQual;
	}

	public float getMinMapq() {
		return minMapq;
	}

	public float getMinVaf() {
		return minVaf;
	}
}
