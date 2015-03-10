/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Breaks contigs down to sub-contigs of length 2*read_length.
 * Filters chunks that do not vary from reference
 * Drops chunks that are wholly contained within another chunk
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ContigChopper {
	
	private int readLength;
	private CompareToReference2 c2r;

	private static int chunkIdx = 0;
	

	//TODO: Merge regions < readLength apart?  Shouldn't matter from accuracy standpoint.
	private Map<String, SAMRecord> chop(SAMFileReader reader, Feature region) {
		// Map of bases to SAMRecords.  Used to dedup reads.
		Map<String, SAMRecord> chunks = new HashMap<String, SAMRecord>();
		
		// Map of positions to base strings.  Used to count noise in a chunk region.
		Map<String, Set<String>> posCounts = new HashMap<String, Set<String>>();

		CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
		
		while (iter.hasNext()) {
			SAMRecord contig = iter.next();
			
			List<SAMRecord> contigChunks = chunkRead(contig);
			for (SAMRecord chunk : contigChunks) {
				chunks.put(chunk.getReadString(), chunk);
			}
		}
		
		iter.close();
		
		return chunks;
	}
	
	public void chopClopDrop(List<Feature> regions, String input, String output) {
		
		SAMFileReader reader = new SAMFileReader(new File(input));
		
		SAMFileHeader header = reader.getFileHeader();
		header.setSortOrder(SortOrder.unsorted);
		
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				header, true, new File(output));
		
		for (Feature region : regions) {
			Map<String, SAMRecord> chunks = chop(reader, region);
			chunks = clop(chunks);
			chunks = drop(chunks);
			
			for (SAMRecord chunk : chunks.values()) {
				out.addAlignment(chunk);
			}
		}
		
		reader.close();
		out.close();
	}
	
	// Filter reads that match reference
	private Map<String, SAMRecord> clop(Map<String, SAMRecord> chunks) {
	
		Map<String, SAMRecord> filtered = new HashMap<String, SAMRecord>();
		
		for (SAMRecord chunk : chunks.values()) {
			int numMismatches = c2r.numMismatches(chunk); 
			
			if (chunk.getCigarLength() != 1) {
				filtered.put(chunk.getReadString(), chunk);
			} else if (chunk.getCigar().getCigarElement(0).getOperator() != CigarOperator.M) {
				filtered.put(chunk.getReadString(), chunk);
			} else {
				if (numMismatches > 0) { 
					filtered.put(chunk.getReadString(), chunk);
				}
			}
		}
		
		return filtered;
	}
	
	// Drop reads that are wholly contained within other reads
	private Map<String, SAMRecord> drop(Map<String, SAMRecord> chunks) {
		
		Set<String> reads = new HashSet<String>();
		
		reads.addAll(chunks.keySet());
		
		for (String read : reads) {
			
			for (String read2 : chunks.keySet()) {
				
				if (read2.length() > read.length()) {
					if (read2.contains(read)) {
						chunks.remove(read);
						break;
					}
				}
			}
		}
		
		return chunks;
	}
		
	private List<SAMRecord> chunkRead(SAMRecord contig) {
		List<SAMRecord> chunks = new ArrayList<SAMRecord>();
		
		SAMRecordUtils.removeSoftClips(contig);
		
		int pos = contig.getAlignmentStart();
		
		int idx = 0;
		
//		pos = (contig.getAlignmentStart() / readLength) * readLength + 1;
		pos = contig.getAlignmentStart();
		
		int lastPos = contig.getAlignmentStart() + contig.getReadLength();
		
		int indelOffset = 0;
		
		int endIdx = -1;
		
//		while (pos < lastPos) {
		while (endIdx < contig.getReadLength()-1) {
			// i.e. 913 / 100 = 9.  9 * 100 + 1 = 901
			
			int actualPos = pos;
			
			if (pos > lastPos - 2*readLength) {
				actualPos = lastPos - 2*readLength;
				pos = lastPos;
			}
			
			actualPos = Math.max(actualPos, contig.getAlignmentStart());
			
			// 0 based
			int startIdx = actualPos - contig.getAlignmentStart();
			endIdx = Math.min(startIdx + readLength*2, contig.getReadLength()-1);
//			int endIdx = Math.min(end - contig.getAlignmentStart(), startIdx+contig.getReadLength()-1); // inclusive
			String bases = contig.getReadString().substring(startIdx, endIdx+1);
			Cigar cigar = SAMRecordUtils.subset(contig.getCigar(), startIdx, endIdx);
			
			SAMRecord chunk = cloneRead(contig);
			// Using a static variable because contigs may overlap multiple regions and cause dupes
			chunk.setReadName(chunk.getReadName() + "_" + chunkIdx++);
			chunk.setAlignmentStart(actualPos + calcIndelOffset(startIdx, contig.getCigar()));
			chunk.setCigar(cigar);
			chunk.setReadString(bases);
			chunk.clearAttributes();
//			chunk.setAttribute("ZZ", contig.getCigarString());
			
			chunk.setAttribute("ZZ", contig.getReferenceName() + ":" + contig.getAlignmentStart() +
					":" + contig.getCigarString()); 
			
			chunks.add(chunk);
			
			pos += readLength;
		}
		
		return chunks;
	}
	
	private int calcIndelOffset(int length, Cigar cigar) {
		int offset = 0;
		int elemLength = 0;
		
		for (CigarElement elem : cigar.getCigarElements()) {
			if (elemLength > length) {
				break;
			} else if ((elemLength == length) && (elem.getOperator() != CigarOperator.D)) {
				break;
			} else if (elemLength + elem.getLength() > length) {
				if (elem.getOperator() == CigarOperator.D) {
					//offset += (length-elemLength);
					offset += elem.getLength();  // Next chunk should be after deletion.  So, add the whole block
				} else if (elem.getOperator() == CigarOperator.I) {
					offset -= (length-elemLength);
				}
				
				if (elem.getOperator() != CigarOperator.D) {
					elemLength += elem.getLength();
				}
			} else {
				if (elem.getOperator() == CigarOperator.D) {
					offset += elem.getLength();
				} else if (elem.getOperator() == CigarOperator.I) {
					offset -= elem.getLength();
				}
				
				if (elem.getOperator() != CigarOperator.D) {
					elemLength += elem.getLength();
				}
			}
		}
		
		return offset;
	}
	
	public void setReadLength(int len) {
		this.readLength = len;
	}
	
	public void setC2R(CompareToReference2 c2r) {
		this.c2r = c2r;
	}
	
	private SAMRecord cloneRead(SAMRecord read) {
		try {
			return (SAMRecord) read.clone();
		} catch (CloneNotSupportedException e) {
			// Infamous "this should never happen" comment here.
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	public static void main(String[] args) throws Exception  {
		CompareToReference2 c2r = new CompareToReference2();
//		c2r.init("/home/lmose/reference/chr6/chr6.fa");
		c2r.init("/home/lmose/reference/chr1/chr1.fa");
		
		RegionLoader loader = new RegionLoader();
//		List<Feature> regions = loader.load("/home/lmose/dev/abra_wxs/3/sim83.gtf");
//		List<Feature> regions = loader.load("/home/lmose/dev/ayc/regions/clinseq5/setd2.gtf");
		List<Feature> regions = loader.load("/home/lmose/dev/ayc/regions/clinseq5/chop3.gtf", false);
//		List<Feature> regions = loader.load("/home/lmose/dev/ayc/regions/clinseq5/uncseq5.gtf");
		
//		regions = ReAligner.splitRegions(regions);
		
		ContigChopper chopper = new ContigChopper();
		chopper.setReadLength(100);
		chopper.setC2R(c2r);

//		chopper.chopClopDrop(regions, "/home/lmose/dev/abra_wxs/chop3/funk.bam",
//				"/home/lmose/dev/abra_wxs/chop3/output3.bam");

		chopper.chopClopDrop(regions, "/home/lmose/dev/abra_wxs/chop3/all_contigs_chim_sorted.bam",
				"/home/lmose/dev/abra_wxs/chop3/output2.bam");

		
//		chopper.chopClopDrop(regions, "/home/lmose/dev/abra_wxs/chop2/fuzzy.bam",
//				"/home/lmose/dev/abra_wxs/chop2/output.bam");
		
//		chopper.chopClopDrop(regions, "/home/lmose/dev/abra_wxs/chop2/t.bam",
//				"/home/lmose/dev/abra_wxs/chop2/to.bam");
		
//		chopper.chopClopDrop(regions, "/home/lmose/dev/abra_wxs/chop2/t2.bam",
//		"/home/lmose/dev/abra_wxs/chop2/to2.bam");


	}
	
	/*
	public static void main(String[] args) throws Exception  {
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("/home/lmose/reference/chr1/1.fa");
		SAMFileReader reader = new SAMFileReader(new File("/home/lmose/dev/abra_wxs/2/chop.bam"));
		
		ContigChopper chopper = new ContigChopper();
		chopper.setReadLength(76);
		chopper.setC2R(c2r);
		
		Feature region = new Feature("1", 110882018, 110884218);
		
		Collection<SAMRecord> chunks = chopper.chop(reader, region);
		chunks = chopper.filterReference(chunks);
		
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				reader.getFileHeader(), true, new File("/home/lmose/dev/abra_wxs/2/out.bam"));
		
		for (SAMRecord chunk : chunks) {
			out.addAlignment(chunk);
		}
		
		out.close();
		reader.close();
	}
	*/
}
