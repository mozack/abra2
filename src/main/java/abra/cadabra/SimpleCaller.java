package abra.cadabra;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

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
	
	private FishersExactTest fishers = new FishersExactTest();
	
	// Indices into orientation counts
	int aIdx = 0;
	int cIdx = 2;
	int tIdx = 4;
	int gIdx = 6;
	int indelIdx = 8;

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
			
			//TODO: Encapsulate in class
			int[] orientationCounts = new int[] { 0,0,0,0,0,0,0,0,0,0 };
			
			if (!c2r.containsChromosome(reads.getChromosome())) {
				System.err.println("Chromosome: [" + reads.getChromosome() + "] not in reference.  Assuming we've reached unaligned pile and stopping.");
				break;
			}
			
			char ref = c2r.containsChromosome(reads.getChromosome()) ? Character.toUpperCase(c2r.getSequence(reads.getChromosome(), reads.getPosition(), 1).charAt(0)) : 'N';
			
			if (ref != 'N') {
			
//				if (reads.getPosition() == 3193842) {
//					System.out.println("yo.");
//				}
				
				int numIndels = 0;
				
				for (SAMRecord read : reads.getReads()) {
					
					if (read.getMappingQuality() >= this.minMapq && !read.getReadUnmappedFlag()) {
					
						BaseInfo baseInfo = getBaseAtReferencePosition(read, reads.getPosition());
						int baseCountIdx = -1;
						
						switch(baseInfo.base) {
							case 'A':
								a++;
								if (baseInfo.isNearEdge) {
									aAtEdge++;
								}
								baseCountIdx = aIdx;
								break;
							case 'C':
								c++;
								if (baseInfo.isNearEdge) {
									cAtEdge++;
								}
								baseCountIdx = cIdx;
								break;
							case 'T':
								t++;
								if (baseInfo.isNearEdge) {
									tAtEdge++;
								}
								baseCountIdx = tIdx;
								break;
							case 'G':
								g++;
								if (baseInfo.isNearEdge) {
									gAtEdge++;
								}
								baseCountIdx = gIdx;
								break;
							default:
								n++;
								break;
						}
						
						if (baseInfo.isIndel) {
							numIndels++;
							baseCountIdx = indelIdx;
						}
						
						if (baseCountIdx > -1) {
							if (!read.getReadNegativeStrandFlag()) {
								// Update forward orientation count
								orientationCounts[baseCountIdx] += 1;
							} else {
								// Update reverse orientation count
								orientationCounts[baseCountIdx+1] += 1;
							}
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
				
				CallInfo callInfo = getAltBaseAndCounts(bases, counts, edgeCounts, ref, reads.getReads().size(), orientationCounts);
				
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
		System.out.println("##FORMAT=<ID=RF,Number=1,Type=String,Description=\"Reference forward strand observations\">");
		System.out.println("##FORMAT=<ID=RR,Number=1,Type=String,Description=\"Reference reverse strand observations\">");
		System.out.println("##FORMAT=<ID=AF,Number=1,Type=String,Description=\"Alternate forward strand observations\">");
		System.out.println("##FORMAT=<ID=AR,Number=1,Type=String,Description=\"Alternate reverse strand observations\">");
		System.out.println("##FORMAT=<ID=FO,Number=1,Type=String,Description=\"Fisher's Exact Test evaluating ref/alt forward/reverse obs\">");
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
		call.append("\t.\t.\tFOO=BAR;\tDP1:DP2:RO:AO:AF1:AF2:RF:RR:AF:AR:FO\t");
		
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
			call.append(':');
			call.append(callInfo.refForward);
			call.append(':');
			call.append(callInfo.refReverse);
			call.append(':');
			call.append(callInfo.altForward);
			call.append(':');
			call.append(callInfo.altReverse);
			call.append(':');
			
			double fs = fishers.twoTailedTest(callInfo.refForward, callInfo.refReverse, callInfo.altForward, callInfo.altReverse);
			double phredFs = -10 * Math.log10(fs);
			call.append(df.format(phredFs));
			
			System.out.println(call.toString());
		}
	}
	
	private CallInfo getAltBaseAndCounts(char[] bases, int[] counts, int[] edgeCounts, char ref, int totalDepth, int[] orientationCounts) {
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
		
		int refIdx = -1;
		int altIdx = -1;
				
		switch (ref) {
			case 'A':
				refIdx = aIdx;
				break;
			case 'C':
				refIdx = cIdx;
				break;
			case 'T':
				refIdx = tIdx;
				break;
			case 'G':
				refIdx = gIdx;
				break;
		}
		
		switch (alt) {
			case 'A':
				altIdx = aIdx;
				break;
			case 'C':
				altIdx = cIdx;
				break;
			case 'T':
				altIdx = tIdx;
				break;
			case 'G':
				altIdx = gIdx;
				break;
		}
		
		int refF = 0;
		int refR = 0;
		int altF = 0;
		int altR = 0;
		
		if (refIdx > -1 && altIdx > -1) {
			refF = orientationCounts[refIdx];
			refR = orientationCounts[refIdx+1];
			altF = orientationCounts[altIdx];
			altR = orientationCounts[altIdx+1];
		}
			
		return new CallInfo(ref, refCount, alt, altCount, altEdgeCount, totalDepth, refF, refR, altF, altR);
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
		int refForward;
		int refReverse;
		int altForward;
		int altReverse;
		
		CallInfo(char ref, int refCount, char alt, int altCount, int altEdgeCount, int totalDepth, int refForward, int refReverse, int altForward, int altReverse) {
			this.ref = ref;
			this.refCount = refCount;
			this.alt = alt;
			this.altCount = altCount;
			this.altEdgeCount = altEdgeCount;
			this.totalDepth = totalDepth;
			this.refForward = refForward;
			this.refReverse = refReverse;
			this.altForward = altForward;
			this.altReverse = altReverse;
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
	
	static class BaseCounts {
		int forward;
		int reverse;
	}
	
	public static void main(String[] args) throws Exception {
		SimpleCaller c = new SimpleCaller();
		
		String reference = args[0];
		String bam = args[1];
		float minAllelicFraction = Float.parseFloat(args[2]);
		int minAltObs = Integer.parseInt(args[3]);
		int minMapq = Integer.parseInt(args[4]);
		int minDistanceFromIndel = Integer.parseInt(args[5]);
	
		c.call(reference, bam, minAllelicFraction, minAltObs, minMapq, minDistanceFromIndel);
		
//		c.call("/home/lmose/reference/chr20/chr20.fa", "/home/lmose/dev/efseq/piotr_test1/calling/k101.sscs.chr20.bam", .003F, 2, 40, 50);
//		c.call("/home/lmose/reference/chr21/chr21.fa", "/home/lmose/dev/efseq/piotr_test1/calling/tiny21.bam", .003F, 2, 40, 50);
		
//		c.call("/home/lmose/reference/chr20/chr20.fa", "/home/lmose/dev/efseq/piotr_test1/calling/k101.chr20.bam", .003F, 2, 40, 50);
		
	}
}
