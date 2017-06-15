package abra.cadabra;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import abra.CompareToReference2;
import abra.Feature;
import abra.Logger;
import abra.SAMRecordUtils;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class CadabraProcessor {

	private static final int MIN_SUPPORTING_READS = 2;
	private static final double MIN_ALLELE_FRACTION = 0.10;
	private static final int MIN_MAPQ = 20;
	
	private Cadabra cadabra;
	private String normalBam;
	private String tumorBam;
	private ReadLocusReader normal;
	private ReadLocusReader tumor;
	private CompareToReference2 c2r;
	private Feature region;
	private int lastPos = 0;
	
	List<SampleCall> sampleRecords = new ArrayList<SampleCall>();
	List<SomaticCall> somaticCalls = new ArrayList<SomaticCall>();
	Map<Integer, SampleCall> normalCalls = new HashMap<Integer, SampleCall>();
	
	CadabraProcessor(Cadabra cadabra, String tumorBam, CompareToReference2 c2r) {
		this.cadabra = cadabra;
		this.tumorBam = tumorBam;
		this.c2r = c2r;
	}
	
	CadabraProcessor(Cadabra cadabra, String normalBam, String tumorBam, CompareToReference2 c2r) {
		this(cadabra, tumorBam, c2r);
		this.normalBam = normalBam;
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
				call.getVaf() >= MIN_ALLELE_FRACTION;
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
					SampleCall normalCall = processLocus(normalReads, true);
					SampleCall tumorCall = processLocus(tumorReads, true);
					
					if (tumorCall.alt != null && tumorCall.alt != Allele.UNK && tumorCall.alleleCounts.get(tumorCall.alt).getCount() >= MIN_SUPPORTING_READS) {
						
						if (normalCall.getVaf()/tumorCall.getVaf() < .2) {
							
							int chromosomeLength = c2r.getChromosomeLength(tumorCall.chromosome);
							String refSeq = "N";
							if (tumorCall.position > 10 && tumorCall.position < chromosomeLength-10) {
								refSeq = c2r.getSequence(tumorCall.chromosome, tumorCall.position-9, 20);
							}
							
							SomaticCall somaticCall = new SomaticCall(normalCall, tumorCall, refSeq);
							somaticCalls.add(somaticCall);
						}
					}
					
					if (normalCall.alt != null && (normalCall.alt.getType() == Allele.Type.DEL || normalCall.alt.getType() == Allele.Type.INS)) {
						normalCalls.put(normalCall.position, normalCall);
					}

					normalReads = normalIter.next();
					tumorReads = tumorIter.next();
				}
				
				if ((count % 1000000) == 0) {
					System.err.println("Position: " + normalReads.getChromosome() + ":" + normalReads.getPosition());
				}
				
				count += 1;
			} else {
				normalReads = normalIter.next();
				tumorReads = tumorIter.next();
			}
		}
		
		// Annotate somatic calls that have overlapping normal indels
		for (SomaticCall call : somaticCalls) {
			int pos = call.tumor.position;
			int normalOverlap = 0;
			
			int stop = pos;
			if (call.tumor.alt.getType() == Allele.Type.DEL) {
				stop += call.tumor.alt.getLength()+1;
			} else {
				stop += 1;
			}
			
			for (int i=pos-100; i<=stop; i++) {
				SampleCall normalCall = normalCalls.get(pos);
				if (normalCall != null) {
					int normalStop = normalCall.alt.getType() == Allele.Type.DEL ? 
							normalCall.position + normalCall.alt.getLength() + 1 : normalCall.position + 1;
					if (normalStop >= pos) {
						normalOverlap += 1;
					}
				}
			}
			
			call.overlappingNormalAF = (float) normalOverlap / (float) call.normal.usableDepth;
		}
		
		normalCalls.clear();
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
				
				if (read.getMappingQuality() < MIN_MAPQ) {
					if (read.getMappingQuality() == 0) {
						tumorMapq0 += 1;
					}
					continue;
				}
				
				if (read.getStringAttribute("YA") == null) {
					// Cap # mismatches in read that can be counted as reference
					// This is done because realigner caps # of mismatches for remapped indel reads.
					// This is needed to remove ref bias
					int editDist = SAMRecordUtils.getEditDistance(read, null);
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
		
		int usableDepth = AlleleCounts.sum(alleleCounts.values());
		
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
			
//			if (altCounts.getCount() >= MIN_SUPPORTING_READS && af >= MIN_ALLELE_FRACTION) {

				double qual = isSomatic ? 0 : calcPhredScaledQuality(refCounts.getCount(), altCounts.getCount(), usableDepth);
				int repeatPeriod = getRepeatPeriod(chromosome, position, alt, altCounts);

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
						usableDepth, qual, repeatPeriod, tumorMapq0, refField, altField, mismatchExceededReads, refSeq);
//			}
		} else {
			String refField = getInsRefField(chromosome, position);
			String altField = ".";
			double qual = 0;
			int rp = 0;
			
			call = new SampleCall(chromosome, position, refAllele, Allele.UNK, alleleCounts, totalDepth, 
					usableDepth, qual, rp, tumorMapq0, refField, altField, mismatchExceededReads, refSeq);
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
		
		public static final String FORMAT = "DP:DP2:AD:AD2:MIRI:MARI:SOR:FS:MQ0:ISPAN:VAF:MER:BB:GT";
		
		String chromosome;
		int position;
		Allele ref;
		Allele alt;
		Map<Allele, AlleleCounts> alleleCounts;
		int totalReads;
		int usableDepth;
		double qual;
		int repeatPeriod;
		int mapq0;
		String refField;
		String altField;
		double fs;
		int mismatchExceededReads;
		HomopolymerRun hrun;
		String context;
		
		SampleCall(String chromosome, int position, Allele ref, Allele alt, Map<Allele, AlleleCounts> alleleCounts, 
				int totalReads, int usableDepth, double qual, int repeatPeriod, int mapq0, String refField, String altField,
				int mismatchExceededReads, String context) {
			this.chromosome = chromosome;
			this.position = position;
			this.ref = ref;
			this.alt = alt;
			this.alleleCounts = alleleCounts;
			this.totalReads = totalReads;
			this.usableDepth = usableDepth;
			this.qual = qual;
			this.repeatPeriod = repeatPeriod;
			this.mapq0 = mapq0;
			this.refField = refField;
			this.altField = altField;
			
			AlleleCounts refCounts = alleleCounts.get(ref);
			AlleleCounts altCounts = alleleCounts.get(alt);
			
			if (refCounts != null && altCounts != null) {
//				this.fs = strandBias(refCounts.getFwd(), refCounts.getRev(), altCounts.getFwd(), altCounts.getRev());
				this.fs = 0;
			}
			
			this.mismatchExceededReads = mismatchExceededReads;
			
			if (context != null) {
				this.hrun = HomopolymerRun.find(context);
				this.context = context;
			}
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
			
			int ispan = altCounts == null ? 0 : altCounts.getMaxReadIdx()-altCounts.getMinReadIdx();
			float vaf = getVaf();
			
			double bbQual = calcPhredScaledQuality(refCounts.getCount(), altCounts.getCount(), usableDepth);
			
			String sampleInfo = String.format("%d:%d:%d,%d:%d,%d:%d:%d:%d,%d,%d,%d:%f:%d:%d:%f:%d:%f:./.", totalReads, usableDepth, refCounts.getCount(), altCounts.getCount(),
					refCounts.getTotalCount(), altCounts.getTotalCount(),
					altCounts.getMinReadIdx(), altCounts.getMaxReadIdx(), refCounts.getFwd(), refCounts.getRev(), altCounts.getFwd(), altCounts.getRev(),
					fs, mapq0, ispan, vaf, mismatchExceededReads, bbQual);

			return sampleInfo;
		}
		
		public String toString() {
			
			//
			// chr1    14397   .       CTGT    C       31.08108108108108       PASS    SOMATIC;CMQ=0;CTX=TAAAAGCACACTGTTGGTTT;REPEAT_PERIOD=1;NNAF=<NNAF>      
			// DP:AD:YM0:YM1:YM:OBS:MIRI:MARI:SOR:MQ0:GT       1092:51,23:0:0:0:23:5:36:0,51,1,22:981:0/1
			String pos = String.valueOf(position);
			String qualStr = String.valueOf(qual);
			
			int hrunLen = hrun != null ? hrun.getLength() : 0;
			char hrunBase = hrun != null ? hrun.getBase() : 'N';
			int hrunPos = hrun != null ? hrun.getPos() : 0;
			
			String info = String.format("REPEAT_PERIOD=%d;HRUN=%d,%c,%d;REF=%s", repeatPeriod,
					hrunLen, hrunBase, hrunPos, context);
			
			String sampleInfo = getSampleInfo(ref, alt);
			
			return String.join("\t", chromosome, pos, ".", refField, altField, qualStr, "PASS", info, SampleCall.FORMAT, sampleInfo);
		}
	}
	
	static double calcFisherExactPhredScaledQuality(int normalRefObs, int normalAltObs, int tumorRefObs, int tumorAltObs) {
		FishersExactTest test = new FishersExactTest();
		// Calc p-value
		double p = test.oneTailedTest(normalRefObs, normalAltObs, tumorRefObs, tumorAltObs);
		
		// Convert to phred scale
		double qual = -10 * Math.log10(p);
		
		// Round to tenths
		qual = (int) (qual * 10);
		qual = qual / 10.0;
		
		return qual;
	}
	
	public static class SomaticCall {
		SampleCall normal;
		SampleCall tumor;
		
		double qual;
		float overlappingNormalAF;
		HomopolymerRun hrun;
		String context;
		
		public SomaticCall(SampleCall normal, SampleCall tumor, String context) {
			this.normal = normal;
			this.tumor = tumor;
			
			int normalRef = normal.alleleCounts.get(tumor.ref) == null ? 0 : normal.alleleCounts.get(tumor.ref).getCount();
			int normalAlt = normal.alleleCounts.get(tumor.alt) == null ? 0 : normal.alleleCounts.get(tumor.alt).getCount();
			
			int tumorRef = tumor.alleleCounts.get(tumor.ref).getCount();
			int tumorAlt = tumor.alleleCounts.get(tumor.alt).getCount();
			
			this.qual = calcFisherExactPhredScaledQuality(normalRef, normalAlt, tumorRef, tumorAlt);
			this.hrun = HomopolymerRun.find(context);
			this.context = context;
		}
		
		public String toString() {
			
			String pos = String.valueOf(tumor.position);
			String qualStr = String.valueOf(qual);
			int hrunLen = hrun != null ? hrun.getLength() : 0;
			char hrunBase = hrun != null ? hrun.getBase() : 'N';
			int hrunPos = hrun != null ? hrun.getPos() : 0;
			
			String info = String.format("REPEAT_PERIOD=%d;ONAF=%f;HRUN=%d,%c,%d;REF=%s", tumor.repeatPeriod, overlappingNormalAF,
					hrunLen, hrunBase, hrunPos, context);
			
			String normalInfo = normal.getSampleInfo(tumor.ref, tumor.alt);
			String tumorInfo = tumor.getSampleInfo(tumor.ref, tumor.alt);
			
			return String.join("\t", tumor.chromosome, pos, ".", tumor.refField, tumor.altField, qualStr, "PASS", info, SampleCall.FORMAT, normalInfo, tumorInfo);
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
	
	private int getRepeatPeriod(String chromosome, int position, Allele indel, AlleleCounts indelCounts) {
		int chromosomeEnd = c2r.getReferenceLength(chromosome);
		int length = Math.min(indel.getLength() * 20, chromosomeEnd-position-2);
		String sequence = c2r.getSequence(chromosome, position+1, length);
		
		String bases;
		if (indel.getType() == Allele.Type.DEL) {
			bases = sequence.substring(0, indel.getLength());
		} else {
			bases = indelCounts.getPreferredInsertBases();
		}
		
		int period = 0;
		
		if (bases.length() > 0) {
			int index = 0;
			while ((index+bases.length() < length) && (bases.equals(sequence.substring(index, index+bases.length())))) {
				period += 1;
				index += bases.length();
			}
		}
		
		return period;
	}
	
	private String getDelRefField(String chromosome, int position, int length) {
		return c2r.getSequence(chromosome, position, length+1);
	}
	
	private String getInsRefField(String chromosome, int position) {
		return c2r.getSequence(chromosome, position, 1);
	}	

	private IndelInfo checkForIndelAtLocus(SAMRecord read, int refPos) {
		IndelInfo elem = null;
		
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
					if (elem.getCigarElement().getOperator() == CigarOperator.I && 
						elem.getCigarElement().getLength() == readElem.getCigarElement().getLength()) {
					
						String insertBases = read.getReadString().substring(readElem.getReadIndex(), readElem.getReadIndex()+readElem.getCigarElement().getLength());
						elem.setInsertBases(insertBases);
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
