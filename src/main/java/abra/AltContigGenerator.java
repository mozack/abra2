package abra;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
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
	
	// Return true if read contains a soft clip element >= 5 bases long with 80% or more of bases exceeding min base qual
	private boolean hasHighQualitySoftClipping(SAMRecord read) {
		
		boolean hasHighQualitySoftClipping = false;
		
		if (read.getCigarLength() > 1) {
			// Check first cigar element
			CigarElement elem = read.getCigar().getCigarElement(0);
			if (elem.getOperator() == CigarOperator.S  && elem.getLength() > minSoftClipLength) {
				int elemStart = 0;
				hasHighQualitySoftClipping = hasHighQualitySoftClipping(read, elemStart, elem.getLength());
			}
			
			// Check last Cigar element
			if (!hasHighQualitySoftClipping) {
				elem = read.getCigar().getCigarElement(read.getCigarLength()-1);
				if (elem.getOperator() == CigarOperator.S  && elem.getLength() > minSoftClipLength) {
					int elemStart = read.getReadLength() - elem.getLength();
					hasHighQualitySoftClipping = hasHighQualitySoftClipping(read, elemStart, elem.getLength());
				}
			}
		}
		
		return hasHighQualitySoftClipping;
	}

	public Collection<String> getAltContigs(List<List<SAMRecordWrapper>> readsList, CompareToReference2 c2r, int readLength) {
		
		List<ScoredContig> softClipContigs = new ArrayList<ScoredContig>();
		
		HashSet<Indel> indels = new HashSet<Indel>();
		
		Map<String, List<ScoredContig>> softClipByPos = new HashMap<String, List<ScoredContig>>();
		
		for (List<SAMRecordWrapper> reads : readsList) {
			for (SAMRecordWrapper readWrapper : reads) {
				SAMRecord read = readWrapper.getSamRecord();
				
				if (read.getMappingQuality() > minMapq) {
					
					if (useObservedIndels) {
						// For now only use indels bracketed by 2 M elements
						List<CigarElement> elems = read.getCigar().getCigarElements();
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
							
							Indel indel = new Indel(type, read.getReferenceName(), read.getAlignmentStart() + elems.get(0).getLength(), elems.get(1).getLength(), insertBases);
							indels.add(indel);
						} else if(SAMRecordUtils.getNumIndels(read) > 1) {
							// Handle read containing multiple indels (create single contig)
							List<Indel> indelComponents = new ArrayList<Indel>();
							List<ReadBlock> readBlocks = SAMRecordUtils.getReadBlocks(read.getCigar(), read.getAlignmentStart());
							
							// Loop through cigar ignoring edge elements
							for (int i=1; i<read.getCigar().getCigarElements().size()-1; i++) {
								CigarElement elem = read.getCigar().getCigarElements().get(i);
								if (elem.getOperator() == CigarOperator.D || elem.getOperator() == CigarOperator.I) {
									ReadBlock block = readBlocks.get(i);
									char type = '0';
									String insertBases = null;
									if (elem.getOperator() == CigarOperator.D) {
										type = 'D';
									} else {
										type = 'I';
										int start = block.getReadPos()-1;
										int stop = start + elem.getLength();
										insertBases = read.getReadString().substring(start, stop);
									}
									
									Indel indel = new Indel(type, read.getReferenceName(), block.getRefPos(), elem.getLength(), insertBases);
									indelComponents.add(indel);
								}
							}
							
							Indel indel = new Indel('C', read.getReferenceName(), indelComponents);
							indels.add(indel);
						}
					}
					
					// Add high quality soft clipped reads
					if (useSoftClippedReads && hasHighQualitySoftClipping(readWrapper.getSamRecord())) {
						
						ScoredContig sc = new ScoredContig((double) SAMRecordUtils.sumBaseQuals(read) / (double) read.getReadLength(), read.getReadString());
						
						if (useConsensusSeq) {
							// Group by position and read length
							int start = readWrapper.getAdjustedAlignmentStart();
							int end = readWrapper.getAdjustedAlignmentEnd();
							String pos = "" + start + ":" + end + ":" + readWrapper.getSamRecord().getReadLength();
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
		List<ScoredContig> filteredContigs = ScoredContig.filter(new ArrayList<ScoredContig>(softClipContigs),  this.maxSoftClipContigs);
		
		Set<String> contigs = new HashSet<String>();
		for (ScoredContig contig : filteredContigs) {
			contigs.add(contig.getContig());
		}
		
		ConsensusSequence cs = new ConsensusSequence();
		for (List<ScoredContig> posList : softClipByPos.values()) {
			contigs.add(cs.buildConsensus(posList));
		}
		
		for (Indel indel : indels) {
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
					seqStart = prev.pos+1;
				}
				String rightSeq = c2r.getSequence(indel.chr, seqStart, readLength);
				seq.append(rightSeq);
				
				contigs.add(seq.toString());
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
		
		Indel(char type, String chr, int pos, int length, String insert) {
			this.type = type;
			this.chr = chr;
			this.pos = pos;
			this.length = length;
			this.insert = insert;
		}
		
		// For complex indels
		Indel(char type, String chr, List<Indel> indels) {
			this.type = type;
			this.chr = chr;
			this.components = indels;
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

	}
}
