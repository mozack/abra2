package abra.cadabra;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import abra.CompareToReference2;
import abra.Feature;
import abra.Logger;
import abra.Pair;
import abra.SAMRecordUtils;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class CadabraProcessor {

	private static final int MIN_SUPPORTING_READS = 2;
	
	private Cadabra cadabra;
	private String normalBam;
	private String tumorBam;
	private ReadLocusReader normal;
	private ReadLocusReader tumor;
	private CompareToReference2 c2r;
	private Feature region;
	private int lastPos = 0;
	CadabraOptions options;
	
	List<SampleCall> sampleRecords = new ArrayList<SampleCall>();
	List<SomaticCall> somaticCalls = new ArrayList<SomaticCall>();
	
	CadabraProcessor(Cadabra cadabra, CadabraOptions options, CompareToReference2 c2r) {
		this.cadabra = cadabra;
		this.options = options;
		this.c2r = c2r;
		this.tumorBam = options.getTumor();
		this.normalBam = options.getNormal();
	}
	
	void process(Feature region) {
		this.region = region;
		this.tumor = new ReadLocusReader(tumorBam, region);
		if (normalBam != null) {
			this.normal = new ReadLocusReader(normalBam, region);
			processSomatic();
		} else {
			processSimple();
		}
	}
		
	private void processSimple() {
		Iterator<ReadsAtLocus> sampleIter = tumor.iterator();
		
		ReadsAtLocus sampleReads = null;
		
		while (sampleIter.hasNext()) {
			
			sampleReads = sampleIter.next();
			SampleCall call = processLocus(sampleReads, false);
			if (call != null && sampleCallExceedsThresholds(call)) {
				sampleRecords.add(call);
			}
		}
		
		this.cadabra.addCalls(region.getSeqname(), sampleRecords);
	}
	
	private boolean sampleCallExceedsThresholds(SampleCall call) {
		return call.alt != null && call.alt != Allele.UNK && call.alleleCounts.get(call.alt).getCount() >= MIN_SUPPORTING_READS &&
				call.getVaf() >= options.getMinVaf() && call.qual >= options.getMinQual();
	}
	
	private int getRepeatLength(int period, String unit, Allele.Type alleleType) {
		if (alleleType == Allele.Type.DEL) {
			// Subtract 1 from deletions as we are looking for reference context
			period = Math.max(period-1, 0);
		}
		
		return period * unit.length();
	}
	
	private void processSomatic() {
		Iterator<ReadsAtLocus> normalIter = normal.iterator();
		Iterator<ReadsAtLocus> tumorIter = tumor.iterator();
		
		ReadsAtLocus normalReads = null;
		ReadsAtLocus tumorReads = null;
		
		int count = 0;
		
		while (normalIter.hasNext() && tumorIter.hasNext()) {
			if (normalReads != null && tumorReads != null) {
				int compare = normalReads.compareLoci(tumorReads, normal.getSamHeader().getSequenceDictionary());
				
				if (compare < 0) {
					normalReads = normalIter.next();
				} else if (compare > 0) {
					tumorReads = tumorIter.next();
				} else {
					SampleCall tumorCall = processLocus(tumorReads, true);
					SampleCall normalCall = processLocus(normalReads, true);
					
					if (tumorCall.alt != null && tumorCall.alt != Allele.UNK && tumorCall.alleleCounts.get(tumorCall.alt).getCount() >= MIN_SUPPORTING_READS) {
						
						int spanEnd = tumorCall.position + getRepeatLength(tumorCall.repeatPeriod, tumorCall.repeatUnit, tumorCall.alt.getType());
						AlleleCounts.setSpanEnd(spanEnd, tumorCall.alleleCounts);
						AlleleCounts.setSpanEnd(spanEnd, normalCall.alleleCounts);
						
						tumorCall.usableDepth = AlleleCounts.sum(tumorCall.alleleCounts.values());
						normalCall.usableDepth = AlleleCounts.sum(normalCall.alleleCounts.values());
						
						if (normalCall.getVaf()/tumorCall.getVaf() < .2) {
							
							int chromosomeLength = c2r.getChromosomeLength(tumorCall.chromosome);
							String refSeq = "N";
							if (tumorCall.position > 10 && tumorCall.position < chromosomeLength-10) {
								refSeq = c2r.getSequence(tumorCall.chromosome, tumorCall.position-9, 20);
							}
							
							SomaticCall somaticCall = new SomaticCall(normalCall, tumorCall, refSeq, options);
							if (somaticCall.qual >= options.getMinQual() && somaticCall.tumor.getVaf() >= options.getMinVaf()) {
								somaticCalls.add(somaticCall);
							}
						}
					}
					
					normalReads = normalIter.next();
					tumorReads = tumorIter.next();
				}
				
				if ((count % 1000000) == 0) {
					Logger.info("Position: " + normalReads.getChromosome() + ":" + normalReads.getPosition());
				}
				
				count += 1;
			} else {
				normalReads = normalIter.next();
				tumorReads = tumorIter.next();
			}
		}
		
		this.cadabra.addSomaticCalls(region.getSeqname(), somaticCalls);
	}
	
	private Character getBaseAtPosition(SAMRecord read, int refPos) {
		int readPos = 0;
		int refPosInRead = read.getAlignmentStart();
		int cigarElementIdx = 0;
		
		while (refPosInRead <= refPos && cigarElementIdx < read.getCigar().numCigarElements() && readPos < read.getReadLength()) {
			CigarElement elem = read.getCigar().getCigarElement(cigarElementIdx++);
			
			switch(elem.getOperator()) {
				case H: //NOOP
					break;
				case S:
				case I:
					readPos += elem.getLength();
					break;
				case D:
				case N:
					refPosInRead += elem.getLength();
					break;
				case M:
					if (refPos < (refPosInRead + elem.getLength())) {
						readPos += refPos - refPosInRead;
						if (readPos < read.getReadLength()) {
							// Found the base.  Return it
							return read.getReadString().charAt(readPos);
						}
					} else {
						readPos += elem.getLength();
						refPosInRead += elem.getLength();
					}
					break;
				default:
					throw new IllegalArgumentException("Invalid Cigar Operator: " + elem.getOperator() + " for read: " + read.getSAMString());					
			}
		}
		
		return null;
	}
	
	private char getRefBase(String chr, int pos) {
		return c2r.getSequence(chr, pos, 1).charAt(0);
	}
	
	
	private Allele getAltIndelAllele(Allele ref, Map<Allele, AlleleCounts> alleleCounts) {
		int maxAlt = 0;
		Allele alt = null;
		
		for (Allele allele : alleleCounts.keySet()) {
			if (allele != ref) {
				AlleleCounts ac = alleleCounts.get(allele);
				if (ac.getCount() > maxAlt && (allele.getType() == Allele.Type.DEL || allele.getType() == Allele.Type.INS)) {
					maxAlt = ac.getCount();
					alt = allele;
				}
			}
		}
		
		return alt;
	}
	
	private SampleCall processLocus(ReadsAtLocus reads, boolean isSomatic) {
		
		SampleCall call = null;
		
		String chromosome = reads.getChromosome();
		int position = reads.getPosition();
				
		if (position > lastPos + 5000000) {
			Logger.info("Processing: %s:%d", chromosome, position);
			lastPos = position;
		}
		
		int tumorMapq0 = 0;
		int mismatchExceededReads = 0;
		int totalDepth = 0;
				
		Map<Allele, AlleleCounts> alleleCounts = new HashMap<Allele, AlleleCounts>();
		
		// Always include ref allele
		char refBase = getRefBase(chromosome, position);
		Allele refAllele = Allele.getAllele(refBase); 
		alleleCounts.put(refAllele, new AlleleCounts());
		
		for (SAMRecord read : reads.getReads()) {
			
			if (!read.getDuplicateReadFlag() && !read.getReadUnmappedFlag() &&
					(read.getFlags() & 0x900) == 0) {

				totalDepth += 1;
				
				if (read.getMappingQuality() < options.getMinMapq()) {
					tumorMapq0 += 1;
					continue;
				}
				
				if (read.getStringAttribute("YA") == null) {
					// Cap # mismatches in read that can be counted as reference
					// This is done because realigner caps # of mismatches for remapped indel reads.
					// This is needed to remove ref bias
					int editDist = SAMRecordUtils.getEditDistance(read, null, false);
					int indelBases = SAMRecordUtils.getNumIndelBases(read);
					int numMismatches = editDist - indelBases;
					
					float mismatchRate = (float) .05;
					if (numMismatches > SAMRecordUtils.getMappedLength(read) * mismatchRate) {
						// Skip this read
						mismatchExceededReads += 1;
						continue;
					}
				}
			
				IndelInfo readElement = checkForIndelAtLocus(read, position);
				
				Allele allele = Allele.UNK;
				
				if (readElement != null) {
					if (readElement.getCigarElement().getOperator() == CigarOperator.D) {
						allele = new Allele(Allele.Type.DEL, readElement.getCigarElement().getLength());
					} else if (readElement.getCigarElement().getOperator() == CigarOperator.I) {
						allele = new Allele(Allele.Type.INS, readElement.getCigarElement().getLength());
					}
				} else {
					Character base     = getBaseAtPosition(read, position);
					Character nextBase = getBaseAtPosition(read, position+1);
					IndelInfo readIndel = checkForIndelAtLocus(read.getAlignmentStart(),
							read.getCigar(), position);
					
					if (readIndel == null && base != null && nextBase != null) {
						allele = Allele.getAllele(base);
					}
				}
				
				if (allele != Allele.UNK) {
					if (!alleleCounts.containsKey(allele)) {
						alleleCounts.put(allele, new AlleleCounts());
					}
					
					AlleleCounts ac = alleleCounts.get(allele);
					
					ac.incrementCount(read);
										
					if (readElement != null) {
						ac.updateReadIdx(readElement.getReadIndex());
					}
					
					if (allele.getType() == Allele.Type.INS) {
						ac.updateInsertBases(readElement.getInsertBases());
					}
				}
			}
		}
		
		// Allow readId sets to be garbage collected.
		for (AlleleCounts counts : alleleCounts.values()) {
			counts.clearReadIds();
		}
		
		Allele alt = getAltIndelAllele(Allele.getAllele(refBase), alleleCounts);
		
		
		
		String refSeq = null;
		if (!isSomatic) {
			int chromosomeLength = c2r.getChromosomeLength(chromosome);
			refSeq = "N";
			if (position > 10 && position < chromosomeLength-10) {
				refSeq = c2r.getSequence(chromosome, position-9, 20);
			}
		}
		
		if (alt != null && (alt.getType() == Allele.Type.DEL || alt.getType() == Allele.Type.INS) && refAllele != Allele.UNK) {
			AlleleCounts altCounts = alleleCounts.get(alt);
			AlleleCounts refCounts = alleleCounts.get(refAllele);
			
			Pair<Integer, String> repeat = getRepeatPeriod(chromosome, position, alt, altCounts);
			
			double qual = 0;
			int usableDepth = 0;
			if (!isSomatic) {
				int repeatLength = getRepeatLength(repeat.getFirst(), repeat.getSecond(), alt.getType());
				AlleleCounts.setSpanEnd(position+repeatLength, alleleCounts);
				usableDepth = AlleleCounts.sum(alleleCounts.values());
				qual = calcPhredScaledQuality(refCounts.getCount(), altCounts.getCount(), usableDepth);
			}

			String refField = "";
			String altField = "";
			if (alt.getType() == Allele.Type.DEL) {
				refField = getDelRefField(chromosome, position, alt.getLength());
				altField = refField.substring(0, 1);
			} else if (alt.getType() == Allele.Type.INS) {
				refField = getInsRefField(chromosome, position);
				altField = refField + getPreferredInsertBases(alt, altCounts);
			}
			
			call = new SampleCall(chromosome, position, refAllele, alt, alleleCounts, totalDepth, 
					usableDepth, qual, repeat.getFirst(), repeat.getSecond(), tumorMapq0, refField, altField, mismatchExceededReads, refSeq, options);
		} else {
			String refField = getInsRefField(chromosome, position);
			String altField = ".";
			double qual = 0;
			int rp = 0;
			String ru = "";
			
			call = new SampleCall(chromosome, position, refAllele, Allele.UNK, alleleCounts, totalDepth, 
					0, qual, rp, ru, tumorMapq0, refField, altField, mismatchExceededReads, refSeq, options);
			
			if (!isSomatic) {
				// Adjust qual score for PCR slippage
				if (options.getStrpThreshold() > 0 && call.repeatPeriod >= options.getStrpThreshold()) {
					// Penalize short tandem repeat expansion / contraction
					qual -= options.getPcrPenalty();
				} else if (options.getHrunThreshold() > 0 && call.hrun != null && call.hrun.getLength() >= options.getHrunThreshold() && Math.abs(call.ref.getLength() - call.alt.getLength())<10) {
					// Filter short indels near homopolymer runs
					qual -= options.getPcrPenalty();
				}
			}
		}
		
		return call;
	}
	
	private String getPreferredInsertBases(Allele allele, AlleleCounts counts) {
		String bases = null;
		if (counts.getPreferredInsertBases().isEmpty()) {
			StringBuffer buf = new StringBuffer();
			for (int i=0; i<allele.getLength(); i++) {
				buf.append('N');
			}
			bases = buf.toString();
		} else {
			bases = counts.getPreferredInsertBases();
		}
		return bases;
	}
	
	public static class SampleCall {
		
		public static final String FORMAT = "GT:DP:DP2:AD:AD2:ROR:LMQ:ISPAN:VAF:MER:FROR";
		
		String chromosome;
		int position;
		Allele ref;
		Allele alt;
		Map<Allele, AlleleCounts> alleleCounts;
		int totalReads;
		int usableDepth;
		double qual;
		int repeatPeriod;
		String repeatUnit;
		int mapq0;
		String refField;
		String altField;
		int mismatchExceededReads;
		HomopolymerRun hrun;
		String context;
		int ispan;
		double fs;
		CadabraOptions options;
		
		SampleCall(String chromosome, int position, Allele ref, Allele alt, Map<Allele, AlleleCounts> alleleCounts, 
				int totalReads, int usableDepth, double qual, int repeatPeriod, String repeatUnit, int mapq0, String refField, String altField,
				int mismatchExceededReads, String context, CadabraOptions options) {
			this.chromosome = chromosome;
			this.position = position;
			this.ref = ref;
			this.alt = alt;
			this.alleleCounts = alleleCounts;
			this.totalReads = totalReads;
			this.usableDepth = usableDepth;
			this.qual = qual;
			this.repeatPeriod = repeatPeriod;
			this.repeatUnit = repeatUnit;
			this.mapq0 = mapq0;
			this.refField = refField;
			this.altField = altField;
			
			AlleleCounts altCounts = alleleCounts.get(alt);
			
			this.mismatchExceededReads = mismatchExceededReads;
			
			if (context != null) {
				this.hrun = HomopolymerRun.find(context);
				this.context = context;
			}
			
			ispan = altCounts == null ? 0 : altCounts.getMaxReadIdx()-altCounts.getMinReadIdx();
			this.options = options;
		}
		
		public float getVaf() {
			float vaf = 0;
			AlleleCounts altCounts = alleleCounts.get(alt);
			if (altCounts != null) {
				vaf = (float) altCounts.getCount() / (float) usableDepth;
			}
			
			return vaf;
		}
		
		public String getSampleInfo(Allele ref, Allele alt) {
			AlleleCounts refCounts = alleleCounts.get(ref);
			AlleleCounts altCounts = alleleCounts.get(alt);
			
			if (refCounts == null) {
				refCounts = AlleleCounts.EMPTY_COUNTS;
			}
			
			if (altCounts == null) {
				altCounts = AlleleCounts.EMPTY_COUNTS;
			}
			
			float vaf = getVaf();
			
			// Calculate phred scaled probability of read orientations occurring by chance
			int refFwd = refCounts.getFwd();
			int refRev = refCounts.getRev();
			int altFwd = altCounts.getFwd();
			int altRev = altCounts.getRev();
			
			FishersExactTest test = new FishersExactTest();
			double fsP = test.twoTailedTest(refFwd, refRev, altFwd, altRev);
			// Use abs to get rid of -0
			this.fs = Math.abs(-10 * Math.log10(fsP));
			
			String sampleInfo = String.format("0/1:%d:%d:%d,%d:%d,%d:%d,%d,%d,%d:%d:%d:%.2f:%d:%.2f", usableDepth, totalReads, 
					refCounts.getCount(), altCounts.getCount(),
					refCounts.getTotalCount(), altCounts.getTotalCount(),
					refCounts.getFwd(), refCounts.getRev(), altCounts.getFwd(), altCounts.getRev(),
					mapq0, ispan, vaf, mismatchExceededReads, fs);

			return sampleInfo;
		}
		
		public String toString() {
			
			String pos = String.valueOf(position);
			String qualStr = String.format("%.2f", qual);
			
			int hrunLen = hrun != null ? hrun.getLength() : 0;
			char hrunBase = hrun != null ? hrun.getBase() : 'N';
			int hrunPos = hrun != null ? hrun.getPos() : 0;
			
			String info = String.format("RP=%d;RU=%s;HRUN=%d,%d;CTX=%s", repeatPeriod, repeatUnit,
					hrunLen, hrunPos, context);
			
			String sampleInfo = getSampleInfo(ref, alt);
			
			String filter = CadabraProcessor.applyFilters(this, null, options, hrunLen, qual);
						
			return String.join("\t", chromosome, pos, ".", refField, altField, qualStr, filter, info, SampleCall.FORMAT, sampleInfo);
		}
	}
	
	static String applyFilters(SampleCall tumor, SampleCall normal, CadabraOptions options, int hrunLen, double qual) {
		String filter = "";
		
		// Filter variants that do not appear in sufficiently varying read positions
		if (tumor.ispan < options.getIspanFilter()) {
			filter += "ISPAN;";
		}
		
		// Too many low mapq reads
		if ((float)tumor.mapq0 / (float)tumor.totalReads > options.getLowMQFilter()) {
			filter += "LOW_MAPQ;";
		} else if (normal != null && (float)normal.mapq0 / (float)normal.totalReads > options.getLowMQFilter()) {
			filter += "LOW_MAPQ;";
		}
		
		if (tumor.fs > options.getFsFilter()) {
			filter += "FS;";
		}
		
		if (qual < options.getQualFilter()) {
			filter += "LOW_QUAL;";
		}
		
		if (filter.equals("")) {
			filter = "PASS";
		}
		
		return filter;
	}
	
	static double calcFisherExactPhredScaledQuality(int normalRefObs, int normalAltObs, int tumorRefObs, int tumorAltObs) {
		FishersExactTest test = new FishersExactTest();
		// Calc p-value
		double p = test.oneTailedTest(normalRefObs, normalAltObs, tumorRefObs, tumorAltObs);
		
		double qual;
		
		if (p <= 0) {
			// Don't allow division by 0 or rounding to negative value.
			qual = 5000.0;
		} else {
			// Convert to phred scale
			qual = -10 * Math.log10(p);
			
			// Round to tenths
			qual = (int) (qual * 10);
			qual = qual / 10.0;
			
			if (qual > 5000.0) {
				qual = 5000.0;
			}
		}
		
		return qual;
	}
	
	public static class SomaticCall {
		SampleCall normal;
		SampleCall tumor;
		
		double qual;
		double fs;
		HomopolymerRun hrun;
		String context;
		CadabraOptions options;
		
		public SomaticCall(SampleCall normal, SampleCall tumor, String context, CadabraOptions options) {
			this.normal = normal;
			this.tumor = tumor;
			
			int normalRef = normal.alleleCounts.get(tumor.ref) == null ? 0 : normal.alleleCounts.get(tumor.ref).getCount();
			int normalAlt = normal.alleleCounts.get(tumor.alt) == null ? 0 : normal.alleleCounts.get(tumor.alt).getCount();
			
			int tumorRef = tumor.alleleCounts.get(tumor.ref).getCount();
			int tumorAlt = tumor.alleleCounts.get(tumor.alt).getCount();
			
			this.qual = calcFisherExactPhredScaledQuality(normalRef, normalAlt, tumorRef, tumorAlt);
						
			this.hrun = HomopolymerRun.find(context);
			
			// Adjust qual score for PCR slippage
			if (options.getStrpThreshold() > 0 && tumor.repeatPeriod >= options.getStrpThreshold()) {
				// Penalize short tandem repeat expansion / contraction
				qual -= options.getPcrPenalty();
			} else if (options.getHrunThreshold() > 0 && hrun != null && hrun.getLength() >= options.getHrunThreshold() && Math.abs(tumor.ref.getLength() - tumor.alt.getLength())<10) {
				// Filter short indels near homopolymer runs
				qual -= options.getPcrPenalty();
			}
			
			this.context = context;
			this.options = options;
		}
		
		public String toString() {
			
			String pos = String.valueOf(tumor.position);
			String qualStr = String.format("%.2f", qual);
			int hrunLen = hrun != null ? hrun.getLength() : 0;
			char hrunBase = hrun != null ? hrun.getBase() : 'N';
			int hrunPos = hrun != null ? hrun.getPos() : 0;
			
			String info = String.format("RP=%d;RU=%s;HRUN=%d,%d;CTX=%s", tumor.repeatPeriod, tumor.repeatUnit,
					hrunLen, hrunPos, context);
			
			String normalInfo = normal.getSampleInfo(tumor.ref, tumor.alt);
			String tumorInfo = tumor.getSampleInfo(tumor.ref, tumor.alt);
			
			String filter = CadabraProcessor.applyFilters(tumor, normal, options, hrunLen, qual);
			
			return String.join("\t", tumor.chromosome, pos, ".", tumor.refField, tumor.altField, qualStr, filter, info, SampleCall.FORMAT, normalInfo, tumorInfo);
		}
	}
	
	static double strandBias(int rf, int rr, int af, int ar) {
		FishersExactTest test = new FishersExactTest();
		double sb = test.twoTailedTest(rf, rf, af, ar);
		return sb;
	}
	
	static double calcPhredScaledQuality(int refObs, int altObs, int dp) {
		return -10 * Math.log10(BetaBinomial.betabinCDF(dp, altObs));
	}
	
	private Pair<Integer, String> getRepeatPeriod(String chromosome, int position, Allele indel, AlleleCounts indelCounts) {
		int chromosomeEnd = c2r.getReferenceLength(chromosome);
		int length = Math.min(indel.getLength() * 100, chromosomeEnd-position-2);
		String sequence = c2r.getSequence(chromosome, position+1, length);
		
		String bases;
		if (indel.getType() == Allele.Type.DEL) {
			bases = sequence.substring(0, indel.getLength());
		} else {
			bases = indelCounts.getPreferredInsertBases();
		}
		
		String repeatUnit = RepeatUtils.getRepeatUnit(bases);
		int period = RepeatUtils.getRepeatPeriod(repeatUnit, sequence);
		
		return new Pair<Integer, String>(period, repeatUnit);
	}
	
	private String getDelRefField(String chromosome, int position, int length) {
		return c2r.getSequence(chromosome, position, length+1);
	}
	
	private String getInsRefField(String chromosome, int position) {
		return c2r.getSequence(chromosome, position, 1);
	}	

	private IndelInfo checkForIndelAtLocus(SAMRecord read, int refPos) {
		IndelInfo elem = null;
		
//		if (refPos == 105243047 && read.getReadName().equals("D7T4KXP1:400:C5F94ACXX:5:2302:20513:30410")) {
//			System.out.println("bar");
//		}
		
		String contigInfo = read.getStringAttribute("YA");
		if (contigInfo != null) {
			// Get assembled contig info.
			String[] fields = contigInfo.split(":");
			int contigPos = Integer.parseInt(fields[1]);
			
			Cigar contigCigar = TextCigarCodec.decode(fields[2]);
			
			// Check to see if contig contains indel at current locus
			elem = checkForIndelAtLocus(contigPos, contigCigar, refPos);
			
			if (elem != null) {
				// Now check to see if this read supports the indel
				IndelInfo readElem = checkForIndelAtLocus(read.getAlignmentStart(),
						read.getCigar(), refPos);
				
				// Allow partially overlapping indels to support contig
				// (Should only matter for inserts)
				if (readElem == null || readElem.getCigarElement().getOperator() != elem.getCigarElement().getOperator()) {
					// Read element doesn't match contig indel
					elem = null;
				} else {
					elem.setReadIndex(readElem.getReadIndex());
					
					// If this read overlaps the entire insert, capture the bases.
					if (elem.getCigarElement().getOperator() == CigarOperator.I) {

						if (elem.getCigarElement().getLength() == readElem.getCigarElement().getLength()) {
					
							String insertBases = read.getReadString().substring(readElem.getReadIndex(), readElem.getReadIndex()+readElem.getCigarElement().getLength());
							elem.setInsertBases(insertBases);
						} else if (readElem.getCigarElement().getLength() < elem.getCigarElement().getLength()) {
							
							int lengthDiff = elem.getCigarElement().getLength() - readElem.getCigarElement().getLength();
							
							if (readElem.getReadIndex() == 0) {
								elem.setReadIndex(readElem.getReadIndex() - lengthDiff);
							} else if (readElem.getReadIndex() == read.getReadLength()-1) {
								elem.setReadIndex(readElem.getReadIndex() + lengthDiff);
							}
						}
					}
				}
			}
		}
		
		return elem;
	}
	
	
	private IndelInfo checkForIndelAtLocus(int alignmentStart, Cigar cigar, int refPos) {
		
		IndelInfo ret = null;
		
		int readIdx = 0;
		int currRefPos = alignmentStart;
		for (CigarElement element : cigar.getCigarElements()) {
			if (element.getOperator() == CigarOperator.M) {
				readIdx += element.getLength();
				currRefPos += element.getLength();
			} else if (element.getOperator() == CigarOperator.I) {
				if (currRefPos == refPos+1) {
					ret = new IndelInfo(element, readIdx);
					break;
				}
				readIdx += element.getLength();
			} else if (element.getOperator() == CigarOperator.D) {
				if (currRefPos == refPos+1) {
					ret = new IndelInfo(element, readIdx);
					break;
				}				
				currRefPos += element.getLength();
			} else if (element.getOperator() == CigarOperator.S) {
				readIdx += element.getLength();
			} else if (element.getOperator() == CigarOperator.N) {
				currRefPos += element.getLength();
			}
			
			if (currRefPos > refPos+1) {
				break;
			}
		}
		
		return ret;
	}
}
