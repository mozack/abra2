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

	public static List<List<SAMRecordWrapper>> getReads(List<String> inputFiles, Feature region, ReAligner realigner) {
		
		List<List<SAMRecordWrapper>> sampleReads = new ArrayList<List<SAMRecordWrapper>>();
		

		for (String input : inputFiles) {
			Set<String> readIds = new HashSet<String>();
			List<SAMRecordWrapper> reads = new ArrayList<SAMRecordWrapper>();
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
				
				boolean shouldAssemble = false;
				SAMRecord read = iter.next();
										
				// If running in paired end mode, drop single end reads
				if (!realigner.isFiltered(read)) {
					// Determine if this read should be included in assembly
					if ( (!read.getDuplicateReadFlag()) && 
						 (!read.getReadFailsVendorQualityCheckFlag()) &&
						 (read.getMappingQuality() >= realigner.getMinMappingQuality() || read.getReadUnmappedFlag()) &&
						 (!readIds.contains(getIdentifier(read)))) {
						
	//					if (read.getReadString().length() > readLength) {
	//						reader.close();
	//						throw new IllegalArgumentException("Maximum read length of: " + readLength +
	//								" exceeded for: " + read.getSAMString());
	//					}
						
						readIds.add(getIdentifier(read));
						shouldAssemble = true;
					} 
					
					reads.add(new SAMRecordWrapper(read, shouldAssemble));
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
