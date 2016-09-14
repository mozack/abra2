package abra;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

public class ReadLoader {

	public static List<List<SAMRecord>> getReads(List<String> inputFiles, Feature region, ReAligner realigner) {
		
		List<List<SAMRecord>> sampleReads = new ArrayList<List<SAMRecord>>();
		

		for (String input : inputFiles) {
			Set<String> readIds = new HashSet<String>();
			List<SAMRecord> reads = new ArrayList<SAMRecord>();
			sampleReads.add(reads);
			
			SAMFileReader reader = new SAMFileReader(new File(input));
			reader.setValidationStringency(ValidationStringency.SILENT);
			
			Iterator<SAMRecord> iter;
			if (region != null) {
				iter = reader.queryOverlapping(region.getSeqname(), (int) region.getStart(), (int) region.getEnd());
			} else {
				iter = reader.iterator();
			}
			
			while (iter.hasNext()) {
				
				SAMRecord read = iter.next();
										
				// Don't allow same read to be counted twice.
				if ( (!realigner.isFiltered(read)) && 
					 (!read.getDuplicateReadFlag()) && 
					 (!read.getReadFailsVendorQualityCheckFlag()) &&
					 (read.getMappingQuality() >= realigner.getMinMappingQuality() || read.getReadUnmappedFlag()) &&
					 (!readIds.contains(getIdentifier(read)))) {
					
//					if (read.getReadString().length() > readLength) {
//						reader.close();
//						throw new IllegalArgumentException("Maximum read length of: " + readLength +
//								" exceeded for: " + read.getSAMString());
//					}
					
					readIds.add(getIdentifier(read));
											
					reads.add(read);
				}
			}
						
			reader.close();
		}
		
		return sampleReads;
	}

	private static String getIdentifier(SAMRecord read) {
		String id = read.getReadName();
		
		if (read.getReadPairedFlag() && read.getSecondOfPairFlag()) {
			id += "_2";
		}
		
		return id;
	}

	
}
