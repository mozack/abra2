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
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class GermlineProcessor {

	private static final int MIN_SUPPORTING_READS = 2;
	private static final int MIN_DISTANCE_FROM_READ_END = 3;
	private static final double MIN_ALLELE_FRACTION = 0.10;
	private static final int MIN_MAPQ = 20;
	
	private Germline cadabra;
	private String tumorBam;
	private ReadLocusReader tumor;
	private CompareToReference2 c2r;
	private Feature region;
	private int lastPos = 0;
	
	List<Call> outputRecords = new ArrayList<Call>();
	
	GermlineProcessor(Germline cadabra, String tumorBam, CompareToReference2 c2r) {
		this.cadabra = cadabra;
		this.tumorBam = tumorBam;
		this.c2r = c2r;
	}
	
	void process(Feature region) {
		this.region = region;
		this.tumor = new ReadLocusReader(tumorBam, region);
		
		process();
	}
		
	private void process() {
		Iterator<ReadsAtLocus> tumorIter = tumor.iterator();
		
		ReadsAtLocus tumorReads = null;
		
		while (tumorIter.hasNext()) {
			
			tumorReads = tumorIter.next();
			processLocus(tumorReads);
//			tumorReads = tumorIter.next();
		}
		
		this.cadabra.addCalls(region.getSeqname(), outputRecords);
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
	
	private boolean matchesReference(SAMRecord read, int refPos) {
		boolean isMatch = false;
	
		if (!read.getReadUnmappedFlag()) {
			Character base = getBaseAtPosition(read, refPos);
			if (base != null) {
				String seq = c2r.getSequence(read.getReferenceName(), refPos, 1);
				isMatch = base.charValue() == seq.charAt(0); 
			}
		}
		
		return isMatch;
	}
	
	// Returns most abundant indel in reads
	private IndelInfo checkForIndelInReadsAtPosition(ReadsAtLocus tumorReads, int position) {
		
		Map<IndelInfo, Integer> indelCounts = new HashMap<IndelInfo, Integer>();
		int totalIndelCounts = 0;
		
		for (SAMRecord read : tumorReads.getReads()) {
			if (!read.getDuplicateReadFlag()) {
				IndelInfo indel = checkForIndelAtLocus(read, position);
				
				int currCount = 0;
				if (indelCounts.containsKey(indel)) {
					currCount = indelCounts.get(indel);
				}
				
				indelCounts.put(indel, currCount + 1);
				totalIndelCounts += 1;
			}
		}
		
		
		IndelInfo maxIndel = null;
		int max = 0;
		for (IndelInfo indel : indelCounts.keySet()) {
			
			if (indelCounts.get(indel) > max) {
				
			}
		}
		
		// TODO: Placeholder to compile
		return null;
	}
	
	private Allele getAltAllele(Allele ref, Map<Allele, AlleleCounts> alleleCounts) {
		int maxAlt = 0;
		Allele alt = null;
		for (Allele allele : alleleCounts.keySet()) {
			if (allele != ref) {
				AlleleCounts ac = alleleCounts.get(allele);
				if (ac.getCount() > maxAlt) {
					maxAlt = ac.getCount();
					alt = allele;
				} else if (ac.getCount() == maxAlt) {
					// TODO: Disambiguate ties (especially inserts of same length)
					alt = null;
				}
			}
		}
		
		return alt;
	}
	
	private void processLocus(ReadsAtLocus tumorReads) {
		String chromosome = tumorReads.getChromosome();
		int position = tumorReads.getPosition();
		
		if (position > lastPos + 5000000) {
			Logger.info("Processing: %s:%d", chromosome, position);
			lastPos = position;
		}
		
		int tumorMapq0 = 0;
		
		// Don't double count overlapping reads.
		// TODO: Determine consensus? 
		Set<String> tumorReadIds = new HashSet<String>();
		
		Map<Allele, AlleleCounts> alleleCounts = new HashMap<Allele, AlleleCounts>();
		
		// Always include ref allele
		char refBase = getRefBase(chromosome, position);
		Allele refAllele = Allele.getAllele(refBase); 
		alleleCounts.put(refAllele, new AlleleCounts());
		
		for (SAMRecord read : tumorReads.getReads()) {
			//TODO: Figure out what to do with non-primary alignments.
			//      Need to reconcile with overlapping read check using read name.
			if (!read.getDuplicateReadFlag() && !read.getReadUnmappedFlag() &&
					!tumorReadIds.contains(read.getReadName()) &&
					(read.getFlags() & 0x900) == 0) {
				
				if (read.getMappingQuality() < MIN_MAPQ) {
					if (read.getMappingQuality() == 0) {
						tumorMapq0 += 1;
					}
					continue;
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
					Character base = getBaseAtPosition(read, position);
					if (base != null) {
						allele = Allele.getAllele(base);
					}
				}
				
				if (!alleleCounts.containsKey(allele)) {
					alleleCounts.put(allele, new AlleleCounts());
				}
				
				AlleleCounts ac = alleleCounts.get(allele);
				ac.incrementCount();
				if (read.getReadNegativeStrandFlag()) {
					ac.incrementRev();
				} else {
					ac.incrementFwd();
				}
				
				if (readElement != null) {
					ac.updateReadIdx(readElement.getReadIndex());
				}
				
				if (allele.getType() == Allele.Type.INS) {
					ac.updateInsertBases(readElement.getInsertBases());
				}
								
				tumorReadIds.add(read.getReadName());
			}
		}
		
		
		Allele alt = getAltAllele(Allele.getAllele(refBase), alleleCounts);
		
		if (alt != null && (alt.getType() == Allele.Type.DEL || alt.getType() == Allele.Type.INS) && refAllele != Allele.UNK) {
			AlleleCounts altCounts = alleleCounts.get(alt);
			AlleleCounts refCounts = alleleCounts.get(refAllele);
			double af = (double) altCounts.getCount() / (double) (altCounts.getCount() + refCounts.getCount());
			
			if (altCounts.getCount() >= MIN_SUPPORTING_READS && af >= MIN_ALLELE_FRACTION) {

				double qual = calcPhredScaledQuality(refCounts.getCount(), altCounts.getCount());
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
				
				Call call = new Call(chromosome, position, refAllele, alt, alleleCounts, tumorReadIds.size(), 
						qual, repeatPeriod, tumorMapq0, refField, altField);
				
//				System.err.println(call);
				outputRecords.add(call);
			}
		}
		
		
//		float tumorFraction = (float) tumorCount / (float) tumorReads.getReads().size();
		
//		float tumorFraction = (float) tumorCount / (float) tumorReadIds.size();
//		
//		if (tumorCount >= MIN_SUPPORTING_READS && hasSufficientDistanceFromReadEnd && tumorFraction >= MIN_ALLELE_FRACTION) {
//			String insertBases = null;
//			if (tumorIndel.getOperator() == CigarOperator.I) {
//				insertBases = getInsertBaseConsensus(insertBasesMap, tumorIndel.getLength());
//			}
//			
//			int repeatPeriod = getRepeatPeriod(chromosome, position, tumorIndel, insertBases);
//			
//			double qual = calcPhredScaledQuality(tumorRefCount, tumorCount);
//			
//			String record = generateRecord(chromosome, position, tumorReads, tumorIndel,
//					tumorCount, tumorRefCount, insertBases, maxContigMapq, mismatch0Count, mismatch1Count, totalMismatchCount, minReadIndex, maxReadIndex,
//					repeatPeriod, qual, tumorRefFwd, tumorRefRev, tumorAltFwd, tumorAltRev,
//					tumorMapq0);
//			
//			this.outputRecords.add(record);
//		}
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
	
	public static class Call {
		String chromosome;
		int position;
		Allele ref;
		Allele alt;
		Map<Allele, AlleleCounts> alleleCounts;
		int totalReads;
		double qual;
		int repeatPeriod;
		int mapq0;
		String refField;
		String altField;
		
		Call(String chromosome, int position, Allele ref, Allele alt, Map<Allele, AlleleCounts> alleleCounts, 
				int totalReads, double qual, int repeatPeriod, int mapq0, String refField, String altField) {
			this.chromosome = chromosome;
			this.position = position;
			this.ref = ref;
			this.alt = alt;
			this.alleleCounts = alleleCounts;
			this.totalReads = totalReads;
			this.qual = qual;
			this.repeatPeriod = repeatPeriod;
			this.mapq0 = mapq0;
			this.refField = refField;
			this.altField = altField;
		}
		
		public String toString() {
			
			AlleleCounts refCounts = alleleCounts.get(ref);
			AlleleCounts altCounts = alleleCounts.get(alt);
			//
			// chr1    14397   .       CTGT    C       31.08108108108108       PASS    SOMATIC;CMQ=0;CTX=TAAAAGCACACTGTTGGTTT;REPEAT_PERIOD=1;NNAF=<NNAF>      
			// DP:AD:YM0:YM1:YM:OBS:MIRI:MARI:SOR:MQ0:GT       1092:51,23:0:0:0:23:5:36:0,51,1,22:981:0/1
			String pos = String.valueOf(position);
			String qualStr = String.valueOf(qual);
			String info = String.format("REPEAT_PERIOD=%d;", repeatPeriod);
			String format = "DP:AD:MIRI:MARI:SOR:MQ0:GT";
			String sample = String.format("%d:%d,%d:%d:%d:%d,%d,%d,%d:%d:0/1", totalReads, refCounts.getCount(), altCounts.getCount(),
					altCounts.getMinReadIdx(), altCounts.getMaxReadIdx(), refCounts.getFwd(), refCounts.getRev(), altCounts.getFwd(), altCounts.getRev(),
					mapq0);
			return String.join("\t", chromosome, pos, ".", refField, altField, qualStr, "PASS", info, format, sample);
		}
	}
	
	// TODO : Not Phred...
	static double calcPhredScaledQuality(int refObs, int altObs) {
		double qual = (double) altObs / ((double) refObs + (double) altObs);
		
		return qual * 100;
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

	
	private int getRepeatPeriod(String chromosome, int position, CigarElement indel, String insertBases) {
		int chromosomeEnd = c2r.getReferenceLength(chromosome);
		int length = Math.min(indel.getLength() * 20, chromosomeEnd-position-2);
		String sequence = c2r.getSequence(chromosome, position+1, length);
		
		String bases;
		if (indel.getOperator() == CigarOperator.D) {
			bases = sequence.substring(0, indel.getLength());
		} else {
			bases = insertBases;
		}
		
		int period = 0;
		int index = 0;
		while ((index+bases.length() < length) && (bases.equals(sequence.substring(index, index+bases.length())))) {
			period += 1;
			index += bases.length();
		}
		
		return period;
	}
	
	private void updateInsertBases(Map<String, Integer> insertBases, String bases) {
		if (insertBases.containsKey(bases)) {
			insertBases.put(bases, insertBases.get(bases) + 1);
		} else {
			insertBases.put(bases, 1);
		}
	}
	
	private String getInsertBaseConsensus(Map<String, Integer> insertBases, int length) {
		int maxCount = -1;
		String maxBases = null;
		
		for (String bases : insertBases.keySet()) {
			int count = insertBases.get(bases);
			if (count > maxCount) {
				maxCount = count;
				maxBases = bases;
			}
		}
		
		if (maxBases == null) {
			StringBuffer buf = new StringBuffer(length);
			for (int i=0; i<length; i++) {
				buf.append('N');
			}
			maxBases = buf.toString();
		}
		
		return maxBases;
	}
	
	private boolean sufficientDistanceFromReadEnd(SAMRecord read, int readIdx) {
		boolean ret = false;
		
		if (readIdx >= MIN_DISTANCE_FROM_READ_END &&
			readIdx <= read.getReadLength()-MIN_DISTANCE_FROM_READ_END-1) {
				ret = true;
		}
		
		return ret;
	}
	
	private String getDelRefField(String chromosome, int position, int length) {
		return c2r.getSequence(chromosome, position, length+1);
	}
	
	private String getInsRefField(String chromosome, int position) {
		return c2r.getSequence(chromosome, position, 1);
	}
	
	private String generateRecord(String chromosome, int position,
			ReadsAtLocus tumorReads, CigarElement indel,
			int tumorObs, int tumorRefObs, String insertBases, int maxContigMapq, int ym0, int ym1, int totalYm,
			int minReadIndex, int maxReadIndex, int repeatPeriod,
			double qual, int tumorRefFwd, int tumorRefRev, int tumorAltFwd, int tumorAltRev,
			int tumorMapq0) {
		
		int tumorDepth = tumorReads.getReads().size();
		
		String context = c2r.getSequence(chromosome, position-10, 20);
		
		StringBuffer buf = new StringBuffer();
		buf.append(chromosome);
		buf.append('\t');
		buf.append(position);
		buf.append("\t.\t");
		
		String ref = ".";
		String alt = ".";
		if (indel.getOperator() == CigarOperator.D) {
			ref = getDelRefField(chromosome, position, indel.getLength());
			alt = ref.substring(0, 1);
		} else if (indel.getOperator() == CigarOperator.I) {
			ref = getInsRefField(chromosome, position);
			alt = ref + insertBases;
		}
		
		buf.append(ref);
		buf.append('\t');
		buf.append(alt);
		buf.append("\t");
		buf.append(qual);
		buf.append("\tPASS\t");
		// <NNAF> is placeholder here
		buf.append("SOMATIC;CMQ=" + maxContigMapq + ";CTX=" + context + ";REPEAT_PERIOD=" + repeatPeriod + ";NNAF=<NNAF>");
		buf.append("\tDP:AD:YM0:YM1:YM:OBS:MIRI:MARI:SOR:MQ0:GT\t");

//		buf.append('\t');
		buf.append(tumorDepth);
		buf.append(':');
		buf.append(tumorRefObs);
		buf.append(',');
		buf.append(tumorObs);
		buf.append(':');
		buf.append(ym0);
		buf.append(':');
		buf.append(ym1);
		buf.append(':');
		buf.append(totalYm);
		buf.append(':');
		buf.append(tumorObs);
		buf.append(':');
		buf.append(minReadIndex);
		buf.append(':');
		buf.append(maxReadIndex);

		buf.append(':');
		buf.append(tumorRefFwd);
		buf.append(',');
		buf.append(tumorRefRev);
		buf.append(',');
		buf.append(tumorAltFwd);
		buf.append(',');
		buf.append(tumorAltRev);
		
		buf.append(':');
		buf.append(tumorMapq0);
		buf.append(":0/1");  // GT placeholder
		
		
		return buf.toString();
	}
	

	private IndelInfo checkForIndelAtLocus(SAMRecord read, int refPos) {
		IndelInfo elem = null;
		
		String contigInfo = read.getStringAttribute("YA");
		if (contigInfo != null) {
			// Get assembled contig info.
			String[] fields = contigInfo.split(":");
			int contigPos = Integer.parseInt(fields[1]);
			
//			Cigar contigCigar = TextCigarCodec.getSingleton().decode(fields[2]);
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
			}
		}
		
		return ret;
	}
	
	private char getReadBase(SAMRecord read, int index) {
		return (char) read.getReadBases()[index];
	}
}
