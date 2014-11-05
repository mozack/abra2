package abra.cadabra;

import java.text.DecimalFormat;
import java.util.Iterator;

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

	public void call(String reference, String bam, float minAltFraction, int minAltObs, int minMapq) throws Exception {
		this.minAltFraction = minAltFraction;
		this.minAltObs = minAltObs;
		this.minMapq = minMapq;
		
		c2r = new CompareToReference2();
		c2r.init(reference);
		
		outputHeader();
		
		this.sample = new ReadLocusReader(bam);
		
		String lastChromosome = "";
		
		Iterator<ReadsAtLocus> sampleIter = sample.iterator();
		
		while (sampleIter.hasNext()) {
			ReadsAtLocus reads = sampleIter.next();
			
			if (!reads.getChromosome().equals(lastChromosome)) {
				System.err.println("Processing chromosome: " + reads.getChromosome());
				lastChromosome = reads.getChromosome();
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
			
			char ref = Character.toUpperCase(c2r.getSequence(reads.getChromosome(), reads.getPosition(), 1).charAt(0));
			
			if (ref != 'N') {
			
//				if (reads.getPosition() == 3193842) {
//					System.out.println("yo.");
//				}
				
				for (SAMRecord read : reads.getReads()) {
					
					if (read.getMappingQuality() >= this.minMapq) {
					
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
					}
				}
				
				char[] bases = { 'A','C','T','G'};
				int[] counts = {a,c,t,g};
				int[] edgeCounts = { aAtEdge, cAtEdge, tAtEdge, gAtEdge };
				
				CallInfo callInfo = getAltBaseAndCounts(bases, counts, edgeCounts, ref, reads.getReads().size());
				
				// Require N number of alt obs not near edge of M block
				if ((callInfo.altCount - callInfo.altEdgeCount) > this.minAltObs && callInfo.altCount > 0) {
					output(reads.getChromosome(), reads.getPosition(), callInfo);
				}
				
//				// Require at least 1 alt obs and # reads indel or clip adjacent < 10% of alt obs
//				if (callInfo.altCount > 1 && callInfo.altCount*.1 > indelOrClipAdjacentCount) {
//					output(reads.getChromosome(), reads.getPosition(), callInfo);
//				}
			}
		}
		
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
	
	private void output(String chromosome, int position, CallInfo callInfo) {
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
					base = (char) read.getReadBases()[readIdx];
					break;
				}
			} else if (element.getOperator() == CigarOperator.I) {
				if (currRefPos == refPos+1) {
					//TODO: Handle insertions
					//ret = new IndelInfo(element, readIdx);
					break;
				}
				readIdx += element.getLength();
			} else if (element.getOperator() == CigarOperator.D) {
				if (currRefPos == refPos+1) {
					//TODO: Handle deletions
					//ret = new IndelInfo(element, readIdx);
					break;
				}				
				currRefPos += element.getLength();
			} else if (element.getOperator() == CigarOperator.S) {
				readIdx += element.getLength();
			}			
		}
		
		return new BaseInfo(Character.toUpperCase(base), isNearEdge);
	}
	
	private boolean isIndelOrClip(CigarElement element) {
		boolean isIndelOrClip = false;
		
		if (element != null) {
			isIndelOrClip = 
					element.getOperator() == CigarOperator.D || element.getOperator() == CigarOperator.I || 
					element.getOperator() == CigarOperator.S || element.getOperator() == CigarOperator.H;
		}
		
		return isIndelOrClip;
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
		
		BaseInfo(char base, boolean isNearEdge) {
			this.base = base;
			this.isNearEdge = isNearEdge;
		}
	}
	
	public static void main(String[] args) throws Exception {
		SimpleCaller c = new SimpleCaller();
		
		String reference = args[0];
		String bam = args[1];
		float minAllelicFraction = Float.parseFloat(args[2]);
		int minAltObs = Integer.parseInt(args[3]);
		int minMapq = Integer.parseInt(args[4]);
	
		c.call(reference, bam, minAllelicFraction, minAltObs, minMapq);
		
//		c.call("/home/lmose/reference/chr20/chr20.fa", "/home/lmose/dev/efseq/piotr_test1/calling/k101.sscs.chr20.bam", 0F, 2);
//		c.call("/home/lmose/reference/chr20/chr20.fa", "/home/lmose/dev/efseq/piotr_test1/calling/k101.chr20.bam", .003F, 2, 40);
		
	}
}
