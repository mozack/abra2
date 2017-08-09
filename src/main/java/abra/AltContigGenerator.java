package abra;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import abra.SAMRecordUtils.ReadBlock;

public class AltContigGenerator {
	
	// TODO: Multiple arbitrary cutoffs here, need to optimize
	// Phred 13 is ~5% error rate
	private int maxSoftClipContigs;
	private int minBaseQual;
	private double softClipFraction;
	private int minSoftClipLength;
	private boolean useObservedIndels;
	private boolean useSoftClippedReads;
	private boolean useConsensusSeq;
	private int minMapq;
	
	public AltContigGenerator(int maxSoftClipContigs, int minBaseQual, int softClipFraction, int minSoftClipLength, 
			boolean useObservedIndels, boolean useSoftClippedReads, boolean useConsensusSeq, int minMapq) {
		this.maxSoftClipContigs = maxSoftClipContigs;
		this.minBaseQual = minBaseQual;
		this.softClipFraction = (double) softClipFraction / 100.0;
		this.minSoftClipLength = minSoftClipLength;
		this.useObservedIndels = useObservedIndels;
		this.useSoftClippedReads = useSoftClippedReads;
		this.useConsensusSeq = useConsensusSeq;
		this.minMapq = minMapq;
	}
	
	private boolean hasHighQualitySoftClipping(SAMRecord read, int start, int length) {
		
		int numHighQualBases = 0;
		int requiredHighQualBases = (int) (softClipFraction * length);
		
		for (int bq = start; bq < start+length; bq++) {
			if (read.getBaseQualities()[bq] >= minBaseQual) {
				numHighQualBases += 1;
				
				if (numHighQualBases >= requiredHighQualBases) {
					return true;
				}
			}
		}
		
		return false;
	}
	
	private boolean hasHighQualitySoftClipping(SAMRecord read, Feature region) {
		
		boolean hasHighQualitySoftClipping = false;
		
		if (read.getCigarLength() > 1) {
			// Check first cigar element
			CigarElement elem = read.getCigar().getCigarElement(0);
			if (elem.getOperator() == CigarOperator.S  && elem.getLength() >= minSoftClipLength && read.getAlignmentStart() >= region.getStart()-read.getReadLength()) {
				int elemStart = 0;
				hasHighQualitySoftClipping = hasHighQualitySoftClipping(read, elemStart, elem.getLength());
			}
			
			// Check last Cigar element
			if (!hasHighQualitySoftClipping) {
				elem = read.getCigar().getCigarElement(read.getCigarLength()-1);
				if (elem.getOperator() == CigarOperator.S  && elem.getLength() >= minSoftClipLength && read.getAlignmentEnd() <= region.getEnd()+read.getReadLength()) {
					int elemStart = read.getReadLength() - elem.getLength();
					hasHighQualitySoftClipping = hasHighQualitySoftClipping(read, elemStart, elem.getLength());
				}
			}
		}
		
		return hasHighQualitySoftClipping;
	}
	
	public Collection<String> getAltContigs(List<List<SAMRecordWrapper>> readsList, CompareToReference2 c2r,
			int readLength, int numJuncPerms, Feature region, List<Variant> knownVariants) {
		
		List<ScoredContig> softClipContigs = new ArrayList<ScoredContig>();
		
		// TODO: Hacky, clean up...
		Map<Indel, Indel> indels = new HashMap<Indel, Indel>();
		List<Indel> indelList = new ArrayList<Indel>();
		
		Map<String, List<ScoredContig>> softClipByPos = new HashMap<String, List<ScoredContig>>();
		
		for (List<SAMRecordWrapper> reads : readsList) {
			
			Set<String> mergedReadIds = new HashSet<String>();
			
			for (SAMRecordWrapper readWrapper : reads) {
				SAMRecord read = readWrapper.getSamRecord();
				
				if (read.getMappingQuality() > minMapq) {
					
					if (useObservedIndels) {
						List<CigarElement> elems = read.getCigar().getCigarElements();
						// Here we require indel to be bracketed by 2 M elements
						if (elems.size() == 3 && 
							elems.get(0).getOperator() == CigarOperator.M && 
							elems.get(2).getOperator() == CigarOperator.M &&
							(elems.get(1).getOperator() == CigarOperator.D || elems.get(1).getOperator() == CigarOperator.I)) {
							
							String insertBases = null;
							char type = '0';
							
							if (elems.get(1).getOperator() == CigarOperator.D) {
								type = 'D';
							} else if (elems.get(1).getOperator() == CigarOperator.I) {
								type = 'I';
								int start = elems.get(0).getLength();
								int stop =  start + elems.get(1).getLength();
								insertBases = read.getReadString().substring(start, stop);
							}
							
							int refStart = read.getAlignmentStart() + elems.get(0).getLength();
							
							// Require indel start to be within region
							if (refStart >= region.getStart() && refStart <= region.getEnd()) {
								Indel indel = new Indel(type, read.getReferenceName(), refStart, elems.get(1).getLength(), insertBases, elems.get(0).getLength(), SAMRecordUtils.sumBaseQuals(read));
								if (indels.containsKey(indel)) {
									indels.get(indel).addReadPosition(elems.get(0).getLength(), SAMRecordUtils.sumBaseQuals(read));
								} else {
									indels.put(indel, indel);
								}
							}
						} else if(SAMRecordUtils.getNumGaps(read) > 1) {
							// Handle read containing multiple indels (create single contig)
							List<Indel> indelComponents = new ArrayList<Indel>();
							List<ReadBlock> readBlocks = SAMRecordUtils.getReadBlocks(read.getCigar(), read.getAlignmentStart());
							int firstIdx = -1;
							
							// Loop through cigar ignoring edge elements
							for (int i=1; i<read.getCigar().getCigarElements().size()-1; i++) {
								CigarElement elem = read.getCigar().getCigarElements().get(i);
								if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.I || elem.getOperator() == CigarOperator.N) {
									ReadBlock block = readBlocks.get(i);
									char type = '0';
									String insertBases = null;
									if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.N) {
										// For purposes of contig sequence generation, we treat introns as deletions here
										type = 'D';
									} else {
										type = 'I';
										int start = block.getReadPos()-1;
										int stop = start + elem.getLength();
										insertBases = read.getReadString().substring(start, stop);
									}
									
									if (firstIdx == -1) {
										firstIdx = block.getReadPos()-1;
									}
																		
									Indel indel = new Indel(type, read.getReferenceName(), block.getRefPos(), elem.getLength(), insertBases, block.getReadPos()-1, SAMRecordUtils.sumBaseQuals(read));
									indelComponents.add(indel);
								}
							}
							
							Indel indel = new Indel('C', read.getReferenceName(), indelComponents, firstIdx, SAMRecordUtils.sumBaseQuals(read));
							if (indels.containsKey(indel)) {
								indels.get(indel).addReadPosition(firstIdx, SAMRecordUtils.sumBaseQuals(read));
							} else {
								indels.put(indel, indel);
							}
						}
					}
					
					// Add high quality soft clipped reads
					if (useSoftClippedReads && readWrapper.shouldAssemble() &&
							hasHighQualitySoftClipping(readWrapper.getSamRecord(), region)) {
						
						// Don't double add merged read sequence
						if (readWrapper.hasMergedSeq()) {
							if (mergedReadIds.contains(read.getReadName())) {
								continue;
							} else {
								mergedReadIds.add(read.getReadName());
							}
						}
						
						ScoredContig sc;
						
						if (readWrapper.hasMergedSeq()) {
							sc = new ScoredContig((double) SAMRecordUtils.sumBaseQuals(readWrapper.getMergedQual()) / (double) readWrapper.getReadLength(), readWrapper.getMergedSeq());
						} else {
							sc = new ScoredContig((double) SAMRecordUtils.sumBaseQuals(read) / (double) read.getReadLength(), read.getReadString());
						}
						
//						ScoredContig sc = new ScoredContig((double) SAMRecordUtils.sumBaseQuals(read) / (double) read.getReadLength(), read.getReadString());
						
						if (useConsensusSeq) {
							// Group by position and read length
							int start = readWrapper.getAdjustedAlignmentStart();
							int end = readWrapper.getAdjustedAlignmentEnd();
							String pos = "" + start + ":" + end + ":" + readWrapper.getReadLength();
							
							if (!softClipByPos.containsKey(pos)) {
								softClipByPos.put(pos, new ArrayList<ScoredContig>());
							}
							softClipByPos.get(pos).add(sc);
							
						} else {
							softClipContigs.add(sc);
						}						
					}
				}
			}
		}

		// Current set of contigs is from soft clipping.  These can grow to be too big.  Downsample if necessary.;
		int maxCombos = 1024;
		int maxSCContigs = numJuncPerms == 0 ? maxSoftClipContigs : Math.min(maxSoftClipContigs, maxCombos/numJuncPerms);
		List<ScoredContig> filteredContigs = ScoredContig.filter(new ArrayList<ScoredContig>(softClipContigs), maxSCContigs);
		if (maxSCContigs != maxSoftClipContigs) {
			Logger.info("MAX_SC_CONTIG\t%s\t%d", region, maxSCContigs);
		}
		
		Set<String> contigs = new HashSet<String>();
		for (ScoredContig contig : filteredContigs) {
			contigs.add(contig.getContig());
		}
		
		ConsensusSequence cs = new ConsensusSequence();
		for (List<ScoredContig> posList : softClipByPos.values()) {
			contigs.add(cs.buildConsensus(posList));
		}
		
		// Do not allow obs indel list to grow too big.
		// Score based upon number of unique read index start positions
		// TODO: parameterize?
		indelList.addAll(indels.values());
		int maxObsIndels = numJuncPerms == 0 ? 8 : Math.min(8, maxCombos/numJuncPerms);
		indelList = Indel.filter(indelList, maxObsIndels);
		
		if (maxObsIndels != 8) {
			Logger.info("MAX_OBS_CONTIG\t%s\t%d", region, maxObsIndels);
		}
		
		for (Indel indel : indelList) {
			if (indel.type == 'D') {
				// Pull in read length sequence from both sides of deletion.
				int leftStart = indel.pos - readLength;
				int rightStart = indel.pos + indel.length;
				String leftSeq = c2r.getSequence(indel.chr, leftStart, readLength);
				String rightSeq = c2r.getSequence(indel.chr, rightStart, readLength);
				String seq = leftSeq + rightSeq;
				
				contigs.add(seq);
			} else if (indel.type == 'I') {
				// Pull in read length sequence on both sides of insertion
				int leftStart = indel.pos - readLength;
				int rightStart = indel.pos;
				String leftSeq = c2r.getSequence(indel.chr, leftStart, readLength);
				String rightSeq = c2r.getSequence(indel.chr, rightStart, readLength);
				String seq = leftSeq + indel.insert + rightSeq;
				
				contigs.add(seq);
			} else if (indel.type == 'C') {
				// Multiple indels in single read
				
				StringBuffer seq = new StringBuffer();
				Indel prev = indel.components.get(0);
				int start = prev.pos - readLength;
				String startSeq = c2r.getSequence(indel.chr, start, readLength);
				seq.append(startSeq);
				if (prev.type == 'I') {
					seq.append(prev.insert);
				}
				
				Indel next = null;

				for (int i=1; i<indel.components.size(); i++) {
					next = indel.components.get(i);
					
					int seqStart = -1;
					if (prev.type == 'D') {
						seqStart = prev.pos + prev.length;
					} else if (prev.type == 'I') {
						seqStart = prev.pos;
					}
					
					String nextSeq = c2r.getSequence(indel.chr, seqStart, next.pos-seqStart);
					seq.append(nextSeq);
					
					if (next.type == 'I') {
						seq.append(next.insert);
					}
					
					prev = next;
				}
				
				// Append readLength's worth of ref sequence to right of last indel
				int seqStart = -1;
				if (prev.type == 'D') {
					seqStart = prev.pos + prev.length;
				} else if (prev.type == 'I') {
					seqStart = prev.pos;
				}
				String rightSeq = c2r.getSequence(indel.chr, seqStart, readLength);
				seq.append(rightSeq);
				
				contigs.add(seq.toString());
			}
		}
		
		// Process known variants
		if (knownVariants != null) {
			for (Variant variant : knownVariants) {
				if (c2r.getReferenceLength(variant.getChr()) > variant.getPosition() + readLength + variant.getRefSpan()) {
					int leftStart = variant.getPosition() - readLength;
					int rightStart = variant.getPosition() + variant.getRef().length();
					String leftSeq = c2r.getSequence(variant.getChr(), leftStart, readLength);
					String varSeq = variant.getAlt();
					String rightSeq = c2r.getSequence(variant.getChr(), rightStart, readLength);
					
					String seq = leftSeq + varSeq + rightSeq;
					contigs.add(seq.toString());
				}
			}
		}
		
		return contigs;
	}
	
	static class Indel {
		char type;
		String chr;
		int pos;
		int length;
		String insert;
		List<Indel> components;
		Set<Integer> readPositions = new HashSet<Integer>();
		int count;
		long baseQualSum;
		
		Indel(char type, String chr, int pos, int length, String insert, int readPos, int readBaseQuals) {
			this.type = type;
			this.chr = chr;
			this.pos = pos;
			this.length = length;
			this.insert = insert;
			readPositions.add(readPos);
			count = 1;
			baseQualSum = readBaseQuals;
		}
		
		// For complex indels
		Indel(char type, String chr, List<Indel> indels, int readPos, int readBaseQuals) {
			this.type = type;
			this.chr = chr;
			this.components = indels;
			readPositions.add(readPos);
			count = 1;
			baseQualSum = readBaseQuals;
		}
		
		void addReadPosition(int readPos, int readBaseQuals) {
			readPositions.add(readPos);
			count += 1;
			
			if (baseQualSum > Long.MAX_VALUE - 10000) {
				// Don't allow overflow...
				baseQualSum = Long.MAX_VALUE;
			} else {
				baseQualSum += readBaseQuals;
			}
		}

		@Override
		public String toString() {
			String str = "";
			if (type == 'C') {
				str += "{";
				for (Indel i : components) {
					str += i + "; ";
				}
				
				str += " : " + this.readPositions.size();
				
				str += "}";
			} else {
				str = String.format("%d%c:%d:%d", this.length, this.type, this.pos, this.readPositions.size());
			}
			
			return str;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((chr == null) ? 0 : chr.hashCode());
			result = prime * result + ((components == null) ? 0 : components.hashCode());
			result = prime * result + ((insert == null) ? 0 : insert.hashCode());
			result = prime * result + length;
			result = prime * result + pos;
			result = prime * result + type;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Indel other = (Indel) obj;
			if (chr == null) {
				if (other.chr != null)
					return false;
			} else if (!chr.equals(other.chr))
				return false;
			if (components == null) {
				if (other.components != null)
					return false;
			} else if (!components.equals(other.components))
				return false;
			if (insert == null) {
				if (other.insert != null)
					return false;
			} else if (!insert.equals(other.insert))
				return false;
			if (length != other.length)
				return false;
			if (pos != other.pos)
				return false;
			if (type != other.type)
				return false;
			return true;
		}

		static class IndelFreqComparator implements Comparator<Indel> {

			@Override
			public int compare(Indel i1, Indel i2) {
				if (i1.readPositions.size() < i2.readPositions.size()) {
					return 1;
				} else if (i1.readPositions.size() > i2.readPositions.size()) {
					return -1;
				} else {
					if (i1.baseQualSum < i2.baseQualSum) {
						return 1;
					} else if (i1.baseQualSum > i2.baseQualSum) {
						return -1;
					} else {
						return 0;
					}
				}
			}
		}
		
		static List<Indel> filter(List<Indel> indels, int maxIndels) {
			if (indels.size() > maxIndels) {
				Logger.debug("Shrinking observed indels from %d to %d", indels.size(), maxIndels);
				Collections.shuffle(indels);
				Collections.sort(indels, new IndelFreqComparator());
				
				// Subset to only the first MAX_CONTIGS
				indels = indels.subList(0, maxIndels);
				
				Logger.debug("Retained: %s", indels);
			}

			return indels;
		}
	}
}
