package abra.cadabra;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

import abra.CompareToReference2;

public class SimpleCaller {
	
	private ReadLocusReader sample;
	private CompareToReference2 c2r;
	private DecimalFormat df = new DecimalFormat("0.000");
	
	private int minAltObs;
	private float minAltFraction;
	private int minMapq;
	private int minDistanceFromIndel;
	
	private List<CachedCall> callCache = new ArrayList<CachedCall>();
	
	private static int MIN_BASE_QUALITY = 20;

	public void call(String reference, String bam, float minAltFraction, int minAltObs, int minMapq, int minDistanceFromIndel) throws Exception {
		this.minAltFraction = minAltFraction;
		this.minAltObs = minAltObs;
		this.minMapq = minMapq;
		this.minDistanceFromIndel = minDistanceFromIndel;
		
		c2r = new CompareToReference2();
		c2r.init(reference);
		
		outputHeader();
		
		this.sample = new ReadLocusReader(bam);
		
		String lastChromosome = "";
		
		int lastIndelPos = -1;
		
		Iterator<ReadsAtLocus> sampleIter = sample.iterator();
		
		while (sampleIter.hasNext()) {
			ReadsAtLocus reads = sampleIter.next();
			
			if (!reads.getChromosome().equals(lastChromosome)) {
				System.err.println("Processing chromosome: " + reads.getChromosome());
				lastChromosome = reads.getChromosome();
				lastIndelPos = -1;
				flushCache();
			}
			
			int a = 0;
			int c = 0;
			int t = 0;
			int g = 0;
			int n = 0;
			
			int aAtEdge = 0;
			int cAtEdge = 0;
			int tAtEdge = 0;
			int gAtEdge = 0;
			
			char ref = c2r.containsChromosome(reads.getChromosome()) ? Character.toUpperCase(c2r.getSequence(reads.getChromosome(), reads.getPosition(), 1).charAt(0)) : 'N';
			
			if (ref != 'N') {
			
//				if (reads.getPosition() == 3193842) {
//					System.out.println("yo.");
//				}
				
				int numIndels = 0;
				
				for (SAMRecord read : reads.getReads()) {
					
					if (read.getMappingQuality() >= this.minMapq && !read.getReadUnmappedFlag()) {
					
						BaseInfo baseInfo = getBaseAtReferencePosition(read, reads.getPosition());
						
						switch(baseInfo.base) {
							case 'A':
								a++;
								if (baseInfo.isNearEdge) {
									aAtEdge++;
								}
								break;
							case 'C':
								c++;
								if (baseInfo.isNearEdge) {
									cAtEdge++;
								}
								break;
							case 'T':
								t++;
								if (baseInfo.isNearEdge) {
									tAtEdge++;
								}
								break;
							case 'G':
								g++;
								if (baseInfo.isNearEdge) {
									gAtEdge++;
								}
								break;
							default:
								n++;
								break;
						}
						
						if (baseInfo.isIndel) {
							numIndels++;
						}
					}
				}
				
				if ((float) numIndels / (float) reads.getReads().size() > this.minAltFraction) {
					// There is indel support at this locus.  Track it so we can filter nearby SNPs.
					lastIndelPos = reads.getPosition();
				}
				
				char[] bases = { 'A','C','T','G'};
				int[] counts = {a,c,t,g};
				int[] edgeCounts = { aAtEdge, cAtEdge, tAtEdge, gAtEdge };
				
				CallInfo callInfo = getAltBaseAndCounts(bases, counts, edgeCounts, ref, reads.getReads().size());
				
				// Require N number of alt obs not near edge of M block
				if ((callInfo.altCount - callInfo.altEdgeCount) > this.minAltObs && callInfo.altCount > 0) {
					output(reads.getChromosome(), reads.getPosition(), callInfo, lastIndelPos);
				}
			}
			
			checkCache(reads.getPosition(), lastIndelPos);
		}
		
		flushCache();
		System.err.println("Done.");
	}
	
	private void outputHeader() {
		System.out.println("##fileformat=VCFv4.1");
		System.out.println("##FORMAT=<ID=DP1,Number=1,Type=String,Description=\"Total read depth\">");
		System.out.println("##FORMAT=<ID=DP2,Number=1,Type=String,Description=\"Ref obs + Alt obs\">");
		System.out.println("##FORMAT=<ID=AO,Number=1,Type=String,Description=\"Alt observations\">");
		System.out.println("##FORMAT=<ID=RO,Number=1,Type=String,Description=\"Ref observations\">");
		System.out.println("##FORMAT=<ID=AF1,Number=1,Type=String,Description=\"Allele Frequency based upon DP1\">");
		System.out.println("##FORMAT=<ID=AF2,Number=1,Type=String,Description=\"Allele Frequency based upon DP2\">");
		System.out.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE");
	}
	
	private void output(String chromosome, int position, CallInfo callInfo, int lastIndelPos) {
		if (position >= lastIndelPos+minDistanceFromIndel) {
			this.callCache.add(new CachedCall(chromosome, position, callInfo));
		}
	}
	
	private void checkCache(int currPosition, int lastIndelPos) {
		Iterator<CachedCall> iter = this.callCache.iterator();
		while (iter.hasNext()) {
			CachedCall call = iter.next();
			if (Math.abs(call.position - lastIndelPos) < minDistanceFromIndel) {
				// Too close to indel.  Discard call.
				iter.remove();
			} else if (currPosition >= call.position + minDistanceFromIndel) {
				// We have progressed minDistanceFromIndel positions away from the call.  Safe to output
				// Output the call and remove from cache
				write(call.chromosome, call.position, call.callInfo);
				iter.remove();
			}
		}
	}
	
	private void flushCache() {
		for (CachedCall call : this.callCache) {
			write(call.chromosome, call.position, call.callInfo);
		}
		
		this.callCache.clear();
	}
		
	private void write(String chromosome, int position, CallInfo callInfo) {
		
		StringBuffer call = new StringBuffer();
		call.append(chromosome);
		call.append('\t');
		call.append(position);
		call.append("\t.\t");
		call.append(callInfo.ref);
		call.append('\t');
		call.append(callInfo.alt);
		call.append("\t.\t.\tFOO=BAR;\tDP1:DP2:RO:AO:AF1:AF2\t");
		
		int depth1 = callInfo.totalDepth;
		int depth2 = callInfo.altCount + callInfo.refCount;
		float vaf1 = (float) callInfo.altCount / (float) callInfo.totalDepth;
		float vaf2 = (float) callInfo.altCount / (float) depth2;
		
		if (vaf2 > minAltFraction) {
			String vaf1Str = df.format(vaf1);
			String vaf2Str = df.format(vaf2);
			
			call.append(depth1);
			call.append(':');
			call.append(depth2);
			call.append(':');
			call.append(callInfo.refCount);
			call.append(':');
			call.append(callInfo.altCount);
			call.append(':');
			call.append(vaf1Str);
			call.append(':');
			call.append(vaf2Str);
			
			System.out.println(call.toString());
		}
	}
	
	private CallInfo getAltBaseAndCounts(char[] bases, int[] counts, int[] edgeCounts, char ref, int totalDepth) {
		int refCount = 0;
		int altCount = 0;
		int altEdgeCount = 0;
		char alt = 'N';
		
		for (int i=0; i<bases.length; i++) {
			if (bases[i] == ref) {
				refCount = counts[i];
			} else {
				if (counts[i] > altCount) {
					alt = bases[i];
					altCount = counts[i];
					altEdgeCount = edgeCounts[i];
				}
			}
		}
		
		return new CallInfo(ref, refCount, alt, altCount, altEdgeCount, totalDepth);
	}
	
	private BaseInfo getBaseAtReferencePosition(SAMRecord read, int refPos) {
		boolean isNearEdge = false;
		boolean isIndel = false;
		int alignmentStart = read.getAlignmentStart();
		Cigar cigar = read.getCigar();
		
		char base = 'N';
		
		int readIdx = 0;
		int currRefPos = alignmentStart;
		
		for (CigarElement element : cigar.getCigarElements()) {
						
			if (element.getOperator() == CigarOperator.M) {
				readIdx += element.getLength();
				currRefPos += element.getLength();
				
				if (currRefPos > refPos) {  // TODO: double check end of read base here...
					
					int offset = currRefPos - refPos;
					
					if ((offset < 3) || (offset+3 >= element.getLength())) {
						// This position is within 3 bases of start/end of alignment or clipping or indel.
						// Tag this base for evaluation downstream.
						isNearEdge = true;
					}
					
					readIdx -= offset;
					if ((readIdx < 0) || (readIdx >= read.getReadBases().length)) {
						System.err.println("Read index out of bounds for read: " + read.getSAMString());
						break;
					}
					
					if (read.getBaseQualities()[readIdx] >= MIN_BASE_QUALITY) {
						base = (char) read.getReadBases()[readIdx];
					}
					break;
				}
			} else if (element.getOperator() == CigarOperator.I) {
				if (currRefPos == refPos) {
					//TODO: Handle insertions
					isIndel = true;
					break;
				}
				readIdx += element.getLength();
			} else if (element.getOperator() == CigarOperator.D) {
				if (refPos >= currRefPos && refPos <= currRefPos+element.getLength()) {
					//TODO: Handle deletions
					isIndel = true;
					break;
				}				
				currRefPos += element.getLength();
			} else if (element.getOperator() == CigarOperator.S) {
				readIdx += element.getLength();
			}			
		}
		
		return new BaseInfo(Character.toUpperCase(base), isNearEdge, isIndel);
	}
	
	static class CallInfo {
		char ref;
		int refCount;
		char alt;
		int altCount;
		int totalDepth;
		int altEdgeCount;
		
		CallInfo(char ref, int refCount, char alt, int altCount, int altEdgeCount, int totalDepth) {
			this.ref = ref;
			this.refCount = refCount;
			this.alt = alt;
			this.altCount = altCount;
			this.altEdgeCount = altEdgeCount;
			this.totalDepth = totalDepth;
		}
	}
	
	static class BaseInfo {
		char base;
		boolean isNearEdge;
		boolean isIndel;
		
		BaseInfo(char base, boolean isNearEdge, boolean isIndel) {
			this.base = base;
			this.isNearEdge = isNearEdge;
			this.isIndel = isIndel;
		}
	}
	
	static class CachedCall {
		String chromosome;
		int position;
		CallInfo callInfo;
		
		CachedCall(String chromosome, int position, CallInfo callInfo) {
			this.chromosome = chromosome;
			this.position = position;
			this.callInfo = callInfo;
		}
	}
	
	public static void main(String[] args) throws Exception {
		SimpleCaller c = new SimpleCaller();
		
		/*
		String reference = args[0];
		String bam = args[1];
		float minAllelicFraction = Float.parseFloat(args[2]);
		int minAltObs = Integer.parseInt(args[3]);
		int minMapq = Integer.parseInt(args[4]);
		int minDistanceFromIndel = Integer.parseInt(args[5]);
	
		c.call(reference, bam, minAllelicFraction, minAltObs, minMapq, minDistanceFromIndel);
		*/
		
//		c.call("/home/lmose/reference/chr20/chr20.fa", "/home/lmose/dev/efseq/piotr_test1/calling/k101.sscs.chr20.bam", .003F, 2, 40, 50);
		c.call("/home/lmose/reference/chr21/chr21.fa", "/home/lmose/dev/efseq/piotr_test1/calling/tiny21.bam", .003F, 2, 40, 50);
		
//		c.call("/home/lmose/reference/chr20/chr20.fa", "/home/lmose/dev/efseq/piotr_test1/calling/k101.chr20.bam", .003F, 2, 40, 50);
		
	}
}
