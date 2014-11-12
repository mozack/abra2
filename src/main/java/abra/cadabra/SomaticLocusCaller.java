package abra.cadabra;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import abra.CompareToReference2;

/**
 * Given a "VCF-like" file of variants to inspect, produces an output file with normal / tumor counts of each variant.
 * SNVs require the alt value to match to be counted.  Indels merely require the existence of an indel at the start postion.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class SomaticLocusCaller {
	
	private List<LocusInfo> loci = new ArrayList<LocusInfo>();
	private CompareToReference2 c2r;

	public void call(String normal, String tumor, String vcf, String reference) throws IOException {
		loadLoci(vcf);
		c2r = new CompareToReference2();
		c2r.init(reference);

		SAMFileReader normalReader = new SAMFileReader(new File(normal));
		normalReader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileReader tumorReader = new SAMFileReader(new File(tumor));
		normalReader.setValidationStringency(ValidationStringency.SILENT);
        
		for (LocusInfo locus : loci) {
			locus.normalCounts = getCounts(normalReader, locus);
			locus.tumorCounts = getCounts(tumorReader, locus);
		}
		
		normalReader.close();
		tumorReader.close();
		
		outputResults();
	}
	
	private void outputResults() {
		System.out.println("##fileformat=VCFv4.1");
		System.out.println("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR");
		
		for (LocusInfo locus : loci) {
			String[] fields = new String[] {
				locus.chromosome,
				String.valueOf(locus.posStart),
				locus.id,
				locus.ref,
				locus.alt,
				"0",
				getFilter(locus),
				"SOMATIC;",
				"DP:RO:AO",
				getCountsStr(locus.normalCounts),
				getCountsStr(locus.tumorCounts)
			};
			
			StringBuffer buf = new StringBuffer();
				
			for (String field : fields) {
				buf.append(field);
				buf.append('\t');
			}

			buf.deleteCharAt(buf.length()-1);
			
			System.out.println(buf.toString());
		}
	}
	
	private String getCountsStr(Counts counts) {
		return "" + counts.depth + ":" + counts.refCount + ":" + counts.altCount;
	}
	
	private String getFilter(LocusInfo locus) {
		String filter = "";
		if (locus.tumorCounts.altCount == 0) {
			filter = "NO_TUMOR_OBS;";
		}
		
		if (locus.normalCounts.altCount >= 0) {
			filter += "NORMAL_OBS";
		}
		
		if (locus.tumorCounts.depth == 0) {
			filter += "NO_TUMOR_COV";
		}
		
		if (locus.normalCounts.depth == 0) {
			filter += "NO_NORMAL_COV";
		}
		
		if (filter.equals("")) {
			filter = "PASS;";
		}
		
		return filter;
	}
	
	private boolean isWithin(int i, int start, int stop) {
		return i >= start && i <= stop;
	}
	
	private boolean hasIndel(SAMRecord read, LocusInfo locus) {
		int readPos = 0;
		int refPosInRead = read.getAlignmentStart();
		int cigarElementIdx = 0;
		
		while (refPosInRead <= locus.posStop && cigarElementIdx < read.getCigar().numCigarElements() && readPos < read.getReadLength()) {
			CigarElement elem = read.getCigar().getCigarElement(cigarElementIdx);
			
			switch(elem.getOperator()) {
				case H: //NOOP
					break;
				case I:
					if (isWithin(refPosInRead, locus.posStart, locus.posStop)) {
						return true;
					}
					// Fall through to next case
				case S:
					readPos += elem.getLength();
					break;
				case D:
					if (isWithin(refPosInRead, locus.posStart, locus.posStop)) {
						return true;
					}
					// Fall through to next case
				case N:
					refPosInRead += elem.getLength();
					break;
				case M:
					readPos += elem.getLength();
					refPosInRead += elem.getLength();
					break;
				default:
					throw new IllegalArgumentException("Invalid Cigar Operator: " + elem.getOperator() + " for read: " + read.getSAMString());					
			}
		}
		
		return false;
	}
	
	private Character getBaseAtPosition(SAMRecord read, int refPos) {
		int readPos = 0;
		int refPosInRead = read.getAlignmentStart();
		int cigarElementIdx = 0;
		
		while (refPosInRead <= refPos && cigarElementIdx < read.getCigar().numCigarElements() && readPos < read.getReadLength()) {
			CigarElement elem = read.getCigar().getCigarElement(cigarElementIdx);
			
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
		
		return 'N';
	}

	private Counts getCounts(SAMFileReader reader, LocusInfo locus) {
		
		int depth = 0;
		int altCount = 0;
		int refCount = 0;
		
		CloseableIterator<SAMRecord> iter = reader.queryOverlapping(locus.chromosome, locus.posStart, locus.posStop);
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			
			Character base = getBaseAtPosition(read, locus.posStart);
			Character refBase = c2r.getSequence(locus.chromosome, locus.posStart, 1).charAt(0);
			
			// Override input with actual reference
			locus.ref = new String(new char[] { refBase });
			
			if (locus.isIndel()) {
				if (hasIndel(read, locus)) {
					altCount += 1;
				} else if (!base.equals('N') && base.equals(refBase)) {
					refCount += 1;
				}
			} else {
				
				if (!base.equals('N') && base.equals(refBase)) {
					refCount += 1;
				} else if (base.equals(locus.alt.charAt(0))) {
					altCount += 1;
				}
			}
		}
		
		iter.close();
		
		return new Counts(refCount, altCount, depth);
	}
	
	private void loadLoci(String vcf) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(vcf));
		
		String line = reader.readLine();
		int count = 0;
		while (line != null) {
			if (!line.startsWith("#")) {
				LocusInfo locus = new LocusInfo(line);
				loci.add(locus);
				count += 1;
			}
			
			line = reader.readLine();
		}
		
		System.err.println("Loaded [" + count + "] loci for inspection");
		
		reader.close();
	}
	
	static class LocusInfo {
		String id;
		String chromosome;
		int posStart;
		int posStop;
		String ref;
		String alt;
		Counts normalCounts;
		Counts tumorCounts;
		
		LocusInfo(String vcfLine) {
			String[] fields = vcfLine.split("\\s+");
			chromosome = fields[0];
			
			id = fields[2];
			ref = fields[3];
			alt = fields[4];
			
			if (fields[1].contains("-")) {
				String[] positions = fields[1].split("-");
				posStart = Integer.parseInt(positions[0]);
				posStop = Integer.parseInt(positions[1]);
			} else {
				posStart = Integer.parseInt(fields[1]);
				posStop = posStart;
				if (isIndel()) {
					posStop += 1;
				}
			}
			
			if (isIndel()) {
				posStart -= 1;
			}
		}
		
		boolean isIndel() {
			return alt.length() != ref.length();
		}
	}
	
	static class Counts {
		int refCount;
		int altCount;
		int depth;
		
		Counts(int refCount, int altCount, int depth) {
			this.refCount = refCount;
			this.altCount = altCount;
			this.depth = depth;
		}
	}
	
	public static void main(String[] args) throws Exception {
		String normal = args[0];
		String tumor = args[1];
		String vcf = args[2];
		String reference = args[3];

		SomaticLocusCaller caller = new SomaticLocusCaller();
		caller.call(normal, tumor, vcf, reference);
	}
}
