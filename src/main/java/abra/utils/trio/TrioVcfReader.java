/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra.utils.trio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import abra.Feature;
import abra.RegionLoader;
import abra.utils.trio.LocusGenotype.Genotype;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;

/**
 * Reader class for SAM or BAM file containing paired reads.
 * Provides the ability to iterate over all mated pairs in the file.
 * Supports multiple mappings for the same read. (i.e. repeating read id's)
 * SAM versus BAM format is determined from the file suffix.  i.e. ".sam" or ".bam"
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class TrioVcfReader implements Iterable<TrioGenotype> {
    
    private LocusGenotype cachedFather;
    private LocusGenotype cachedMother;
    private LocusGenotype cachedChild;
    
    
    private SAMRecord lastRead;
    private int lineCnt = 0;
    
    private BufferedReader father;
    private BufferedReader mother;
    private BufferedReader child;
    
    private Map<String, Long> chromosomeOrder;
    
    static final String NO_CHROMOSOME = "";
    
	private Feature currentRegion;
	private Iterator<Feature> regionIter;
	
	public enum Caller { FREEBAYES, GATK, ISAAC };
	
	private Caller caller;
    
    public TrioVcfReader(String chromosomeFile, String father, String mother, String child, Caller caller) throws IOException {
    	
    	this.caller = caller;
//    	GtfLoader loader = new GtfLoader();
//		List<Feature> regions = loader.load(regionsGtf);
//    	
//        regionIter = regions.iterator();
//        if (regionIter.hasNext()) {
//        	currentRegion = regionIter.next();
//        }
    	
    	initChromosomeOrder(chromosomeFile);
    	
		//GZIPInputStream fatherGzip = new GZIPInputStream(new FileInputStream(father));
		//GZIPInputStream motherGzip = new GZIPInputStream(new FileInputStream(mother));
		//GZIPInputStream childGzip = new GZIPInputStream(new FileInputStream(child));
    	
		//this.father = new BufferedReader(new InputStreamReader(fatherGzip));
//		this.mother = new BufferedReader(new InputStreamReader(motherGzip));
//		this.child = new BufferedReader(new InputStreamReader(childGzip));
		
		
		
//		String l = this.father.readLine();
//		while (l != null) {
//			System.out.println(l);
//			l = this.father.readLine();
//		}
		
    	this.father = new BufferedReader(new FileReader(father));
    	this.mother = new BufferedReader(new FileReader(mother));
    	this.child = new BufferedReader(new FileReader(child));
    }
    
    private void initChromosomeOrder(String chromosomeFile) throws IOException {
    	
    	chromosomeOrder = new HashMap<String, Long>();
    	
    	BufferedReader reader = new BufferedReader(new FileReader(chromosomeFile));

    	long order = 0;
    	String chromosome = reader.readLine();
    	
    	while (chromosome != null) {
    		chromosome = chromosome.trim();
    		chromosomeOrder.put(chromosome, order);
    		order += 1000000000;
    		chromosome = reader.readLine();
    	}
    	
    	chromosomeOrder.put(NO_CHROMOSOME, order);
    	
    	reader.close();
    }
    
    private long getChromosomeVal(String chr) {
    	return chromosomeOrder.get(chr);
    }
    
    public void close() throws IOException {
        father.close();
        mother.close();
        child.close();
    }
    
    private LocusGenotype parseLine(String line) {
    	LocusGenotype lgt = null;
    	
    	if (caller == Caller.FREEBAYES) {
    		lgt = parseFreebayesLine(line);
    	} else if (caller == Caller.GATK) {
    		lgt = parseGatkLine(line);
    	} else if (caller == Caller.ISAAC) {
    		lgt = parseIsaacLine(line);
    	} else {
    		throw new IllegalArgumentException("Invalid caller specified...");
    	}
    	
    	return lgt;
    }
    
    private LocusGenotype parseIsaacLine(String line) {
//    	if (line == null) {
//    		return LocusGenotype.REFERENCE;
//    	}
    	try {
	    	String[] fields = line.split("\\s+");
	    	
	    	String chr = fields[0];
	    	int pos = Integer.parseInt(fields[1]);
	    	int end = pos;
	    	
	    	String ref = fields[3];
	    	String[] alt = fields[4].split(",");
	    	double qual = Double.parseDouble(fields[5]);
//	    	String filter = fields[6];
//	    	String info = fields[7];
	    	
//	    	int endIdx = info.indexOf("END=");
//	    	if (endIdx >= 0) {
//	    		int endStopIdx = info.indexOf(";");
//	    		end = Integer.parseInt(info.substring(endIdx+4, endStopIdx));
//	    	}
	    	
	    	Genotype gt = null;
	    	double altAlleleFreq = 0.0;
	    	
	    	
	    	String gtString = fields[9].split(":")[0];
	    	
	    	if (gtString.equals("0/0")) {
				gt = Genotype.REF_REF;
			} else if (gtString.equals("0/1") || gtString.equals("1/0")) {
				gt = Genotype.REF_ALT1;
			} else if (gtString.equals("1/1")) {
				gt = Genotype.ALT1_ALT1;
			} else if (gtString.equals("1/2")) {
				gt = Genotype.ALT1_ALT2;
			} else {			
				gt = Genotype.UNK;
			}
	    	
	    	if ((gt == Genotype.REF_ALT1) && ((ref.length() > 1 && alt[0].length() == 1) || (ref.length() == 1 && alt[0].length() > 1))) {
	    		//String ad = fields[9].split(":")[1];
	    		int adIdx = -1;
	    		if (fields[8].contains("DPI")) {
	    			adIdx = 4;
	    		} else {
	    			adIdx = 5;
	    		}
	    		String ad = fields[9].split(":")[adIdx];
	    		double ad1 = Integer.parseInt(ad.split(",")[0]);
	    		double ad2 = Integer.parseInt(ad.split(",")[1]);
//	    		double ad2 = Integer.parseInt(fields[9].split(":")[4]);
	    		altAlleleFreq = ad2/(ad1+ad2);
	    		
//	    		if (altAlleleFreq < .2) {
//	    			return LocusGenotype.REFERENCE;
//	    		}
	    		
	    	} else if (gt == Genotype.ALT1_ALT1) {
	    		altAlleleFreq = 1.0;
	    	}
	    	
	    	return new LocusGenotype(chr, pos, end, ref, alt[0], alt.length > 1 ? alt[1] : null, gt, qual, altAlleleFreq);
    	} catch (NumberFormatException e) {
    		e.printStackTrace();
    		throw new RuntimeException("Error processing: " + line);
    	}
    }
    
    private LocusGenotype parseGatkLine(String line) {
//    	if (line == null) {
//    		return LocusGenotype.REFERENCE;
//    	}
    	try {
	    	String[] fields = line.split("\\s+");
	    	
	    	String chr = fields[0];
	    	int pos = Integer.parseInt(fields[1]);
	    	int end = pos;
	    	
	    	String ref = fields[3];
	    	String[] alt = fields[4].split(",");
	    	double qual = Double.parseDouble(fields[5]);
//	    	String filter = fields[6];
//	    	String info = fields[7];
	    	
//	    	int endIdx = info.indexOf("END=");
//	    	if (endIdx >= 0) {
//	    		int endStopIdx = info.indexOf(";");
//	    		end = Integer.parseInt(info.substring(endIdx+4, endStopIdx));
//	    	}
	    	
	    	Genotype gt = null;
	    	double altAlleleFreq = 0.0;
	    	
	    	
	    	String gtString = fields[9].split(":")[0];
	    	
	    	if (gtString.equals("0/0")) {
				gt = Genotype.REF_REF;
			} else if (gtString.equals("0/1")) {
				gt = Genotype.REF_ALT1;
			} else if (gtString.equals("1/1")) {
				gt = Genotype.ALT1_ALT1;
			} else if (gtString.equals("1/2")) {
				gt = Genotype.ALT1_ALT2;
			} else {			
				gt = Genotype.UNK;
			}
	    	
	    	if ((gt == Genotype.REF_ALT1) && ((ref.length() > 1 && alt[0].length() == 1) || (ref.length() == 1 && alt[0].length() > 1))) {
	    		String ad = fields[9].split(":")[1];
	    		double depth = Integer.parseInt(fields[9].split(":")[2]);
	    		double ad1 = Integer.parseInt(ad.split(",")[0]);
	    		double ad2 = Integer.parseInt(ad.split(",")[1]);
//	    		double ad2 = Integer.parseInt(fields[9].split(":")[4]);
	    		altAlleleFreq = ad2/(ad1+ad2);
	    		
//	    		if (altAlleleFreq < .2) {
//	    			return LocusGenotype.REFERENCE;
//	    		}
	    		
	    	} else if (gt == Genotype.ALT1_ALT1) {
	    		altAlleleFreq = 1.0;
	    	}
	    	
	    	return new LocusGenotype(chr, pos, end, ref, alt[0], alt.length > 1 ? alt[1] : null, gt, qual, altAlleleFreq);
    	} catch (NumberFormatException e) {
    		e.printStackTrace();
    		throw new RuntimeException("Error processing: " + line);
    	}
    }

    
    private LocusGenotype parseFreebayesLine(String line) {
//    	if (line == null) {
//    		return LocusGenotype.REFERENCE;
//    	}
    	try {
	    	String[] fields = line.split("\\s+");
	    	
	    	String chr = fields[0];
	    	int pos = Integer.parseInt(fields[1]);
	    	int end = pos;
	    	
	    	String ref = fields[3];
	    	String[] alt = fields[4].split(",");
	    	double qual = Double.parseDouble(fields[5]);
//	    	String filter = fields[6];
//	    	String info = fields[7];
	    	
//	    	int endIdx = info.indexOf("END=");
//	    	if (endIdx >= 0) {
//	    		int endStopIdx = info.indexOf(";");
//	    		end = Integer.parseInt(info.substring(endIdx+4, endStopIdx));
//	    	}
	    	
	    	Genotype gt = null;
	    	double altAlleleFreq = 0.0;
	    	
	    	if (line.contains("APRIM")) {
	    		String gtString = fields[9];
	    		if (gtString.equals("0|0")) {
					gt = Genotype.REF_REF;
				} else if (gtString.equals("0|1")) {
					gt = Genotype.REF_ALT1;
				} else if (gtString.equals("1|0")) {
					gt = Genotype.REF_ALT1;
				} else if (gtString.equals("1|1")) {
					gt = Genotype.ALT1_ALT1;
				} else if (gtString.equals("1|2")) {
					gt = Genotype.ALT1_ALT2;
				} else if (gtString.equals("2|1")) {
					gt = Genotype.ALT1_ALT2;
				} else {			
					gt = Genotype.UNK;
				}
	    		
	    		altAlleleFreq = 1.0;
	    	} else {
	    	
		    	String gtString = fields[9].split(":")[0];
		    	
		    	
		    	
		    	if (gtString.equals("0/0")) {
					gt = Genotype.REF_REF;
				} else if (gtString.equals("0/1") || gtString.equals("1/0")) {
					gt = Genotype.REF_ALT1;
				} else if (gtString.equals("1/1")) {
					gt = Genotype.ALT1_ALT1;
				} else if (gtString.equals("1/2")) {
					gt = Genotype.ALT1_ALT2;
				} else {			
					gt = Genotype.UNK;
				}
		    	
		    	
		    	
	//	    	if ((gt == Genotype.REF_ALT1) && (fields[8].contains("DPI"))) {
		    	if ((gt == Genotype.REF_ALT1) && ((ref.length() > 1 && alt[0].length() == 1) || (ref.length() == 1 && alt[0].length() > 1))) {
//		    		String ad = fields[9].split(":")[4];
//		    		double depth = Integer.parseInt(fields[9].split(":")[1]);
		    		double ad1 = Integer.parseInt(fields[9].split(":")[3]);
		    		double ad2 = Integer.parseInt(fields[9].split(":")[5]);
	//	    		String[] alleleDepth = ad.split(",");
	//	    		double ad1 = Integer.parseInt(alleleDepth[0]);
	//	    		double ad2 = Integer.parseInt(alleleDepth[1]);
		    		altAlleleFreq = ad2/(ad1+ad2);
		    		
//		    		if (altAlleleFreq < .2) {
//		    			return LocusGenotype.REFERENCE;
//		    		}

//		    		altAlleleFreq = ad2/(depth);
		    		
	//	    		if (altAlleleFreq < .2) {
	//	    			gt = Genotype.REF_REF;
	//	    		}
		    	} else if (gt == Genotype.ALT1_ALT1) {
		    		altAlleleFreq = 1.0;
		    	}
	    	}
	    	
	    	return new LocusGenotype(chr, pos, end, ref, alt[0], alt.length > 1 ? alt[1] : null, gt, qual, altAlleleFreq);
    	} catch (NumberFormatException e) {
    		e.printStackTrace();
    		throw new RuntimeException("Error processing: " + line);
    	}
    }
    
    private TrioGenotype getNextTrio()  {
    	LocusGenotype father = null;
    	LocusGenotype mother = null;
    	LocusGenotype child  = null;
    	
    	try {
    		String fatherLine = null;
    		String motherLine = null;
    		String childLine = null;
    		
    		if (cachedFather != null) {
    			father = cachedFather;
    			cachedFather = null;
    		} else {
	    		fatherLine = this.father.readLine();
	    		while (fatherLine != null && fatherLine.startsWith("#")) {
	    			fatherLine = this.father.readLine();
	    		}
    		}
    		
    		if (cachedMother != null) {
	    		mother = cachedMother;
	    		cachedMother = null;
	    	} else {
	    		motherLine = this.mother.readLine();
	    		while(motherLine != null && motherLine.startsWith("#")) {
	    			motherLine = this.mother.readLine();
	    		}
	    	}
    		
    		if (cachedChild != null) {
	    		child = cachedChild;
	    		cachedChild = null;
	    	} else {
	    		childLine = this.child.readLine();
	    		while(childLine != null && childLine.startsWith("#")) {
	    			childLine = this.child.readLine();
	    		}
	    	}
    		
        	if (fatherLine == null && motherLine == null && childLine == null) {
        		// We're at the end of all 3 vcfs.
        		return null;
        	}
    		
        	if (fatherLine != null) {
        		father = parseLine(fatherLine);
        	}

        	if (motherLine != null) {
        		mother = parseLine(motherLine);
        	}
        	
        	if (childLine != null) {
        		child = parseLine(childLine);
        	}

    	} catch (IOException e) {
    		e.printStackTrace();
    		throw new RuntimeException(e);
    	}
    	
//    	if (father.getPos() == 9996638) {
//    		System.out.println("here");
//    	}
    	

//    	if (father.getPos() == 47555134) {
//    		System.out.println("here");
//    	}
    	
    	List<Long> ends = new ArrayList<Long>();
    	if (father != null) {
    		ends.add(getEnd(father));
    	}
    	if (mother != null) {
    		ends.add(getEnd(mother));
    	}
    	if (child != null) {
    		ends.add(getEnd(child));
    	}
    	
    	long minEnd = Long.MAX_VALUE;
    	for (long end : ends) {
    		if (end < minEnd) {
    			minEnd = end;
    		}
    	}
    	
    	if (child != null && getEnd(child) > minEnd) {
    		cachedChild = child;
    		child = null;
    	}
    	
    	if (mother != null && getEnd(mother) > minEnd) {
    		cachedMother = mother;
    		mother = null;
    	}
    	
    	if (father != null && getEnd(father) > minEnd) {
    		cachedFather = father;
    		father = null;
    	}
    	
    	/*
    	String minChromosome = father.getChromosome();
    	if (chromosomeOrder.get(mother.getChromosome()) < chromosomeOrder.get(minChromosome)) {
    		minChromosome = mother.getChromosome();
    	}
    	if (chromosomeOrder.get(child.getChromosome()) < chromosomeOrder.get(minChromosome)) {
    		minChromosome = child.getChromosome();
    	}
    	
    	int minPos = Integer.MAX_VALUE;
    	if (father.getChromosome().equals(minChromosome)) {
    		if (father.getPos() < minPos) {
    			minPos = father.getPos();
    		}
    	}
    	
    	if (mother.getChromosome().equals(minChromosome)) {
    		if (mother.getPos() < minPos) {
    			minPos = mother.getPos();
    		}
    	}
    	
    	if (child.getChromosome().equals(minChromosome)) {
    		if (child.getPos() < minPos) {
    			minPos = child.getPos();
    		}
    	}
    	
    	if (!father.getChromosome().equals(minChromosome) || father.getPos() != minPos) {
    		if (father != LocusGenotype.REFERENCE) {
	    		cachedFather = father;
	    		father = LocusGenotype.REFERENCE;
    		}
    	}
    	
    	if (!mother.getChromosome().equals(minChromosome) || mother.getPos() != minPos) {
    		if (mother != LocusGenotype.REFERENCE) {
	    		cachedMother = mother;
	    		mother = LocusGenotype.REFERENCE;
    		}
    	}
    	
    	if (!child.getChromosome().equals(minChromosome) || child.getPos() != minPos) {
    		if (child != LocusGenotype.REFERENCE) {
	    		cachedChild = child;
	    		child = LocusGenotype.REFERENCE;
    		}
    	}
    	*/
    	
    	if (father == null) {
    		father = LocusGenotype.REFERENCE;
    	}
    	
    	if (mother == null) {
    		mother = LocusGenotype.REFERENCE;
    	}
    	
    	if (child == null) {
    		child = LocusGenotype.REFERENCE;
    	}
    	
    	return new TrioGenotype(father, mother, child);
    }
   
    private long getEnd(LocusGenotype lg) {
    	return getChromosomeVal(lg.getChromosome()) + lg.getEnd();
    }

    public Iterator<TrioGenotype> iterator() {
        return new TrioGenotypeIterator(this);
    }
   
	private boolean isInRegion(SAMFileHeader header, String chromosome, int start, int stop) {
//		boolean isInRegion = currentRegion.overlapsRead(read);
		
		while ((currentRegion != null) &&
			   (!currentRegion.overlaps(chromosome, start, stop)) &&
			   (!isRegionBeyondRead(currentRegion, chromosome, start, stop))) {
			
			if (regionIter.hasNext()) {
				currentRegion = (Feature) regionIter.next();
			} else {
				currentRegion = null;
			}
		}
		
		return (currentRegion != null) && (currentRegion.overlaps(chromosome, start, stop));
	}
	
	private boolean isRegionBeyondRead(Feature region, String chromosome, int start, int stop) {
		boolean isRegionBeyond = false;

		
//		if (regionChrIdx > readChrIdx) {
//			isRegionBeyond = true;
//		} else if (regionChrIdx < readChrIdx) {
//			isRegionBeyond = false;
//		} else {
//			isRegionBeyond = (region.getStart() > read.getAlignmentEnd());
//		}
		
		isRegionBeyond = (region.getStart() > stop);
		
		return isRegionBeyond;
	}
    
    
    private static class TrioGenotypeIterator implements Iterator<TrioGenotype> {

        private TrioVcfReader reader;
        private TrioGenotype nextTrio = null;
        
        TrioGenotypeIterator(TrioVcfReader reader) {
            this.reader = reader;
        }
        
        @Override
        public boolean hasNext() {
            if (nextTrio == null) {
            	nextTrio = reader.getNextTrio();
            }
            
            return nextTrio != null;
        }

        @Override
        public TrioGenotype next() {
        	TrioGenotype trio = nextTrio;
        	nextTrio = null;
            return trio;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove not supported for TrioGenotypeIterator.");
        }
    }

    /*
    public static void main(String[] args) throws Exception {
    	String chromosomes = "/home/lmose/dev/abra/trio/chromosomes.txt";
    	String father = "/home/lmose/dev/abra/trio/father.vcf";
    	String mother = "/home/lmose/dev/abra/trio/mother.vcf";
    	String child = "/home/lmose/dev/abra/trio/child.vcf";
    	
    	
    	TrioVcfReader rdr = new TrioVcfReader(chromosomes, father, mother, child);
    	
    	int l = 0;
    	for (TrioGenotype gt : rdr) {
    		if (gt.hasVariant()) {
    			System.out.println(gt.summary());
    		}
    		
//    		if (l++ > 500) break;
    	}
    }
    */
}
