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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

public class GermlineProcessor {

	private static final int MIN_SUPPORTING_READS = 2;
	private static final int MIN_DISTANCE_FROM_READ_END = 3;
	private static final double MIN_TUMOR_FRACTION = 0.02;
	private static final int MIN_TUMOR_MAPQ = 20;
	
	private Germline cadabra;
	private String tumorBam;
	private ReadLocusReader tumor;
	private CompareToReference2 c2r;
	private Feature region;
	
	List<String> outputRecords = new ArrayList<String>();
	
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
			tumorReads = tumorIter.next();
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
	
	private void processLocus(ReadsAtLocus tumorReads) {
		String chromosome = tumorReads.getChromosome();
		int position = tumorReads.getPosition();
		
		CigarElement tumorIndel = null;
		int tumorCount = 0;
		int tumorRefCount = 0;
		int mismatch0Count = 0;
		int mismatch1Count = 0;
		int totalMismatchCount = 0;
		boolean hasSufficientDistanceFromReadEnd = false;
		int maxContigMapq = 0;
		int minReadIndex = Integer.MAX_VALUE;
		int maxReadIndex = Integer.MIN_VALUE;
		int tumorRefFwd = 0;
		int tumorRefRev = 0;
		int tumorAltFwd = 0;
		int tumorAltRev = 0;
		int tumorMapq0 = 0;
		
		Map<String, Integer> insertBasesMap = new HashMap<String, Integer>();
		
		// Don't double count overlapping reads.
		// TODO: Determine consensus? 
		Set<String> tumorReadIds = new HashSet<String>();
		
		for (SAMRecord read : tumorReads.getReads()) {
			//TODO: Figure out what to do with non-primary alignments.
			//      Need to reconcile with overlapping read check using read name.
			if (!read.getDuplicateReadFlag() && !read.getReadUnmappedFlag() && (read.getFlags() & 0x900) == 0) {
				
				if (read.getMappingQuality() < MIN_TUMOR_MAPQ) {
					if (read.getMappingQuality() == 0) {
						tumorMapq0 += 1;
					}
					continue;
				}
			
				IndelInfo readElement = checkForIndelAtLocus(read, position);
				
				if (readElement != null) {
					Integer ymInt = (Integer) read.getAttribute("YM");
					if (ymInt != null) {
						int ym = ymInt;
						if (ym == 0) {
							mismatch0Count++;
						}
						if (ym <= 1) {
							mismatch1Count++;
						}
						totalMismatchCount += ym;
					}
				} else if (matchesReference(read, position)) {
					// TODO: Check for agreement in overlapping read pairs for tumor.
					if (!tumorReadIds.contains(read.getReadName())) {
						tumorRefCount += 1;
						if (read.getReadNegativeStrandFlag()) {
							tumorRefRev += 1;
						} else {
							tumorRefFwd += 1;
						}
					}
				}
				
				if (tumorIndel == null && readElement != null) {
					tumorIndel = readElement.getCigarElement();
					tumorCount = 1;
					if (read.getReadNegativeStrandFlag()) {
						tumorAltRev += 1;
					} else {
						tumorAltFwd += 1;
					}
//					maxContigMapq = Math.max(maxContigMapq, read.getIntegerAttribute(ReadAdjuster.CONTIG_QUALITY_TAG));
					maxContigMapq = 0;
					if (readElement.getInsertBases() != null) {
						updateInsertBases(insertBasesMap, readElement.getInsertBases());
					}
					minReadIndex = readElement.getReadIndex() < minReadIndex ? readElement.getReadIndex() : minReadIndex;
					maxReadIndex = readElement.getReadIndex() > maxReadIndex ? readElement.getReadIndex() : maxReadIndex;
				} else if (tumorIndel != null && readElement != null) {
					if (tumorIndel.equals(readElement.getCigarElement())) {
						// Increment tumor indel support count
						// Do not allow single fragment to contribute multiple reads
						// TODO: Identify consensus
						if (!tumorReadIds.contains(read.getReadName())) {
							tumorCount += 1;
							if (read.getReadNegativeStrandFlag()) {
								tumorAltRev += 1;
							} else {
								tumorAltFwd += 1;
							}
						}
//						maxContigMapq = Math.max(maxContigMapq, read.getIntegerAttribute(ReadAdjuster.CONTIG_QUALITY_TAG));
						maxContigMapq = 0;
						if (readElement.getInsertBases() != null) {
							updateInsertBases(insertBasesMap, readElement.getInsertBases());
						}
						minReadIndex = readElement.getReadIndex() < minReadIndex ? readElement.getReadIndex() : minReadIndex;
						maxReadIndex = readElement.getReadIndex() > maxReadIndex ? readElement.getReadIndex() : maxReadIndex;
	
					} else {
						// We will not deal with multiple indels at a single locus for now.
						tumorIndel = null;
						tumorCount = 0;
						break;
					}
				}
				
				if (!hasSufficientDistanceFromReadEnd && tumorIndel != null && readElement != null && readElement.getCigarElement().equals(tumorIndel)) {
					hasSufficientDistanceFromReadEnd = sufficientDistanceFromReadEnd(read, readElement.getReadIndex());
				}
				
				tumorReadIds.add(read.getReadName());
			}
		}
		
//		float tumorFraction = (float) tumorCount / (float) tumorReads.getReads().size();
		float tumorFraction = (float) tumorCount / (float) tumorReadIds.size();
		
		if (tumorCount >= MIN_SUPPORTING_READS && hasSufficientDistanceFromReadEnd && tumorFraction >= MIN_TUMOR_FRACTION) {
			String insertBases = null;
			if (tumorIndel.getOperator() == CigarOperator.I) {
				insertBases = getInsertBaseConsensus(insertBasesMap, tumorIndel.getLength());
			}
			
			int repeatPeriod = getRepeatPeriod(chromosome, position, tumorIndel, insertBases);
			
			double qual = calcPhredScaledQuality(tumorRefCount, tumorCount);
			
			String record = generateRecord(chromosome, position, tumorReads, tumorIndel,
					tumorCount, tumorRefCount, insertBases, maxContigMapq, mismatch0Count, mismatch1Count, totalMismatchCount, minReadIndex, maxReadIndex,
					repeatPeriod, qual, tumorRefFwd, tumorRefRev, tumorAltFwd, tumorAltRev,
					tumorMapq0);
			
			this.outputRecords.add(record);
		}
	}
	
	// TODO : Not Phred...
	static double calcPhredScaledQuality(int refObs, int altObs) {
		double qual = (double) altObs / ((double) refObs + (double) altObs);
		
		return qual * 100;
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
		buf.append("\tDP:AD:YM0:YM1:YM:OBS:MIRI:MARI:SOR:MQ0\t");

		buf.append('\t');
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
