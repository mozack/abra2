package abra.cadabra;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import abra.CompareToReference2;
import abra.Feature;
import abra.Pair;
import abra.SAMRecordUtils;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class SimpleAlleleCounter {

	private static final int MIN_MAPQ = 20;
	private static final int MIN_BASEQ = 20;
	
	private String bam;
	private ReadLocusReader sample;
	private CompareToReference2 c2r;
	private String inputVcf;
	private List<InputVariant> inputVariants;
	
	List<SampleCall> sampleRecords = new ArrayList<SampleCall>();
	
	SimpleAlleleCounter(CompareToReference2 c2r, String bam, String vcf) {
		this.c2r = c2r;
		this.bam = bam;
		this.inputVcf = vcf;
	}
	
	void run() throws IOException {
		loadInputVariants();
		for (InputVariant variant : inputVariants) {
			process(variant);
		}
	}
	
	void process(InputVariant variant) throws IOException {
		Feature region = new Feature(variant.getChrom(), variant.getPos(), variant.getPos());
		sample = new ReadLocusReader(bam, region);
		processSimple(variant);
		sample.close();
	}
	
	private void loadInputVariants() throws IOException {
		inputVariants = new ArrayList<InputVariant>();
		
		BufferedReader reader = new BufferedReader(new FileReader(inputVcf));
		String line = reader.readLine();
		
		while (line != null) {
			if (!line.startsWith("#") && !line.trim().equals("")) {
				inputVariants.add(InputVariant.create(line));
			}
			
			line = reader.readLine();
		}
		
		reader.close();
	}
		
	private void processSimple(InputVariant variant) {
		Iterator<ReadsAtLocus> sampleIter = sample.iterator();
		
		ReadsAtLocus sampleReads = null;
		
		SampleCall call = null;
		
		while (sampleIter.hasNext() && call == null) {
			
			sampleReads = sampleIter.next();
			call = processLocus(sampleReads, variant);
		}
		
		if (call == null) {
			// Create empty call
			Allele ref = Allele.getAllele(variant.getRef().charAt(0));
			Allele alt = variant.getAllele();
			
			call = SampleCall.emptyCall(variant.getChrom(), variant.getPos(), ref, alt, variant.getRef(), variant.getAlt());
		}
		
		System.out.println(call);
	}
	
//	private boolean sampleCallExceedsThresholds(SampleCall call) {
//		return call.alt != null && call.alt != Allele.UNK && call.alleleCounts.get(call.alt).getCount() >= MIN_SUPPORTING_READS &&
//				call.getVaf() >= options.getMinVaf() && call.qual >= options.getMinQual();
//	}
	
	private int getRepeatLength(int period, String unit, Allele.Type alleleType) {
		if (alleleType == Allele.Type.DEL) {
			// Subtract 1 from deletions as we are looking for reference context
			period = Math.max(period-1, 0);
		} else if (alleleType != Allele.Type.INS) {
			period = 0;
		} 
		
		return period * unit.length();
	}
	
	// Returns Pair of <base, quality(phred)>
	private Pair<Character,Character> getBaseAtPosition(SAMRecord read, int refPos) {
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
							return new Pair<Character, Character>(read.getReadString().charAt(readPos), read.getBaseQualityString().charAt(readPos));
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
	
	private SampleCall processLocus(ReadsAtLocus reads, InputVariant variant) {
		// Always false here
		boolean isSomatic = false;
		
		SampleCall call = null;
		
		String chromosome = reads.getChromosome();
		int position = reads.getPosition();
		
		// Only process positions 
		if (!variant.getChrom().equals(chromosome) || variant.getPos() != position) {
			return null;
		}
		
		int tumorMapq0 = 0;
		int mismatchExceededReads = 0;
		int totalDepth = 0;
				
		Map<Allele, AlleleCounts> alleleCounts = new HashMap<Allele, AlleleCounts>();
		
		// Always include ref allele
		char refBase = getRefBase(chromosome, position);
		Allele refAllele = Allele.getAllele(refBase); 
		alleleCounts.put(refAllele, new AlleleCounts());
		
		// Add input variant allele
		alleleCounts.put(variant.getAllele(), new AlleleCounts());
		
		for (SAMRecord read : reads.getReads()) {
			
			if (!read.getDuplicateReadFlag() && !read.getReadUnmappedFlag() &&
					(read.getFlags() & 0x900) == 0) {
				
				totalDepth += 1;
				
				if (read.getMappingQuality() < MIN_MAPQ) {
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
					// Pair in format <base, quality>
					Pair<Character, Character> base = getBaseAtPosition(read, position);
					Pair<Character, Character> nextBase = getBaseAtPosition(read, position+1);
					IndelInfo readIndel = checkForIndelAtLocus(read.getAlignmentStart(),
							read.getCigar(), position);
					
					if (readIndel == null && base != null && nextBase != null && base.getSecond()-'!' >= MIN_BASEQ) {
						allele = Allele.getAllele(base.getFirst());
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
		
//		Allele alt = getAltIndelAllele(Allele.getAllele(refBase), alleleCounts);
		Allele alt = variant.getAllele();
		
		String refSeq = null;
		if (!isSomatic) {
			int chromosomeLength = c2r.getChromosomeLength(chromosome);
			refSeq = "N";
			if (position > 10 && position < chromosomeLength-10) {
				refSeq = c2r.getSequence(chromosome, position-9, 20);
			}
		}
		
//		if (alt != null && (alt.getType() == Allele.Type.DEL || alt.getType() == Allele.Type.INS) && refAllele != Allele.UNK) {
		if (alt != null && refAllele != Allele.UNK) {
			AlleleCounts altCounts = alleleCounts.get(alt);
			AlleleCounts refCounts = alleleCounts.get(refAllele);
			
			Pair<Integer, String> repeat = getRepeatPeriod(chromosome, position, alt, variant.getAlt());
			
			double qual = 0;
			int usableDepth = 0;

			int repeatLength = getRepeatLength(repeat.getFirst(), repeat.getSecond(), alt.getType());
			AlleleCounts.setSpanEnd(position+repeatLength, alleleCounts);
			usableDepth = AlleleCounts.sum(alleleCounts.values());
			qual = calcPhredScaledQuality(refCounts.getCount(), altCounts.getCount(), usableDepth);

			String refField = variant.getRef();
			String altField = variant.getAlt();
			
			// TODO: Check preferred insert bases against input variant!!!
			
			String altInsert = null;
			
			if (variant.getAllele().getType() == Allele.Type.INS) {
				altInsert = refField + getPreferredInsertBases(alt, altCounts);
				
			}
			
//			if (alt.getType() == Allele.Type.DEL) {
//				refField = getDelRefField(chromosome, position, alt.getLength());
//				altField = refField.substring(0, 1);
//			} else if (alt.getType() == Allele.Type.INS) {
//				refField = getInsRefField(chromosome, position);
//				altField = refField + getPreferredInsertBases(alt, altCounts);
//			}
			
			call = new SampleCall(chromosome, position, refAllele, alt, alleleCounts, totalDepth, 
					usableDepth, qual, repeat.getFirst(), repeat.getSecond(), tumorMapq0, refField, altField, mismatchExceededReads, refSeq, altInsert);
		} else {
			String refField = getInsRefField(chromosome, position);
			String altField = ".";
			double qual = 0;
			int rp = 0;
			String ru = "";
			
			call = new SampleCall(chromosome, position, refAllele, Allele.UNK, alleleCounts, totalDepth, 
					0, qual, rp, ru, tumorMapq0, refField, altField, mismatchExceededReads, refSeq, "");
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
		String altInsert;
		
		static SampleCall emptyCall(String chromosome, int position, Allele ref, Allele alt, String refField, String altField) {
			SampleCall call = new SampleCall();
			call.chromosome = chromosome;
			call.position = position;
			call.ref = ref;
			call.alt = alt;		
			call.alleleCounts = new HashMap<Allele, AlleleCounts>();
			call.alleleCounts.put(ref, new AlleleCounts());
			call.alleleCounts.put(alt, new AlleleCounts());
			call.refField = refField;
			call.altField = altField;
			
			return call;
		}
		
		private SampleCall() {
		}
		
		SampleCall(String chromosome, int position, Allele ref, Allele alt, Map<Allele, AlleleCounts> alleleCounts, 
				int totalReads, int usableDepth, double qual, int repeatPeriod, String repeatUnit, int mapq0, String refField, String altField,
				int mismatchExceededReads, String context, String altInsert) {
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
			this.altInsert = altInsert;
			
			AlleleCounts altCounts = alleleCounts.get(alt);
			
			this.mismatchExceededReads = mismatchExceededReads;
			
			if (context != null) {
				this.hrun = HomopolymerRun.find(context);
				this.context = context;
			}
			
			ispan = altCounts == null ? 0 : altCounts.getMaxReadIdx()-altCounts.getMinReadIdx();
		}
		
		public float getVaf() {
			float vaf = 0;
			AlleleCounts altCounts = alleleCounts.get(alt);
			if (altCounts != null && usableDepth > 0) {
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
			
			String info;
			if (totalReads == 0) {
				// Skip empty call
				info = ".";
			}
			else if (altInsert != null && altField.length() > 1 && !altInsert.equals(altField) && alleleCounts.get(alt).getCount() > 0) {
				// Record info plus alternative inserted sequence (same length, but base mismatches with input variant)
				info = String.format("RP=%d;RU=%s;HRUN=%d,%d;CTX=%s;ALT_INSERT=%s", repeatPeriod, repeatUnit,
						hrunLen, hrunPos, context, altInsert.substring(1));
			} else {
				info = String.format("RP=%d;RU=%s;HRUN=%d,%d;CTX=%s", repeatPeriod, repeatUnit,
					hrunLen, hrunPos, context);
			}
			
			String sampleInfo = getSampleInfo(ref, alt);
						
			return String.join("\t", chromosome, pos, ".", refField, altField, qualStr, ".", info, SampleCall.FORMAT, sampleInfo);
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
	
	private Pair<Integer, String> getRepeatPeriod(String chromosome, int position, Allele indel, String altString) {
		int chromosomeEnd = c2r.getReferenceLength(chromosome);
		int length = Math.min(indel.getLength() * 100, chromosomeEnd-position-2);
		String sequence = c2r.getSequence(chromosome, position+1, length);
		
		String bases;
		if (indel.getType() == Allele.Type.DEL) {
			bases = sequence.substring(0, indel.getLength());
		} else {
			bases = altString.substring(1);
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
	
	static class InputVariant {
		
		private String chrom;
		private int pos;
		private String ref;
		private String alt;
		private Allele allele;
		
		static InputVariant create(String str) {
			String[] fields = str.split("\\s");
			String chrom = fields[0];
			int pos = Integer.parseInt(fields[1]);
			String ref = fields[3];
			String alt = fields[4];
			int length = 1;
			
			if (ref.length() != 1 && alt.length() != 1) {
				// Only supporting simple variant representations for now.
				throw new UnsupportedOperationException("At least one of the REF and ALT fields must be of length 1.");
			}
			
			Allele allele = Allele.UNK;
			
			if (ref.length() > alt.length()) {
				length = ref.length() - alt.length();
				allele = new Allele(Allele.Type.DEL, length);
			} else if (alt.length() > ref.length()) {
				length = alt.length() - ref.length();
				allele = new Allele(Allele.Type.INS, length);
			} else {
				allele = Allele.getAllele(alt.charAt(0));
			}
			
			return new InputVariant(chrom, pos, ref, alt, allele);			
		}
		
		private InputVariant(String chrom, int pos, String ref, String alt, Allele allele) {
			this.chrom = chrom;
			this.pos = pos;
			this.ref = ref;
			this.alt = alt;
			this.allele = allele;
		}

		public String getChrom() {
			return chrom;
		}

		public int getPos() {
			return pos;
		}

		public String getRef() {
			return ref;
		}

		public String getAlt() {
			return alt;
		}
		
		public Allele getAllele() {
			return allele;
		}
	}
	
	public static void main(String[] args) throws Exception {
		String bam = "/home/lmose/dev/mc3/allele_counter/TCGA-D1-A163/TCGA-D1-A163.star.abra2.mc3.bam";
//		String vcf = "/home/lmose/dev/mc3/allele_counter/TCGA-D1-A163/TCGA-D1-A163.maf.mc3.vcf";
//		String vcf = "/home/lmose/dev/mc3/allele_counter/TCGA-D1-A163/TCGA-D1-A163.maf.mc3.chr1.vcf";
		String vcf = "t5.vcf";
		String ref = "/home/lmose/dev/reference/hg19/1.fa";
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(ref);
		SimpleAlleleCounter sac = new SimpleAlleleCounter(c2r, bam, vcf);
		sac.run();
	}
}
