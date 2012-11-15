package abra;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

/**
 * Iterates over a Sam or Bam file returning a List of SAMRecords for each unique read name.
 * Assumes the input file is sorted by read name.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class SamMultiMappingReader implements Iterable<List<SAMRecord>> {
	
	private SAMFileReader inputSam;
	private Iterator<SAMRecord> iter;
	private SAMRecord cachedRead;
	private int numRecords = 0;

    public SamMultiMappingReader(String filename) {
        File inputFile = new File(filename);
        
        inputSam = new SAMFileReader(inputFile);
        inputSam.setValidationStringency(ValidationStringency.SILENT);
  
        iter = inputSam.iterator();
    }
    
    public SAMFileHeader getFileHeader() {
    	return inputSam.getFileHeader();
    }
    
    private boolean hasMoreReads() {
    	return cachedRead != null || iter.hasNext();
    }
    
    private SAMRecord getNextRead() {
    	SAMRecord read = null;
    	
    	if (cachedRead != null) {
    		read = cachedRead;
    		cachedRead = null;
    	} else {
    		read = iter.next();
    		numRecords++;
    		
    		if ((numRecords % 1000000) == 0) {
    			System.out.println("Processed: " + numRecords + " records.");
    		}
    	}
    	
    	return read;
    }
    
    private List<SAMRecord> getNextReadList() {
    	List<SAMRecord> reads = new ArrayList<SAMRecord>();
    	
    	SAMRecord read = null;
    	String readName = null;
    	
    	if (hasMoreReads()) {
    		read = getNextRead();
    		readName = read.getReadName();
    		
    		reads.add(read);
    		
    		while (hasMoreReads() && read.getReadName().equals(readName)) {
    			read = getNextRead();
    			
    			if (read.getReadName().equals(readName)) {
    				reads.add(read);
    			}
    		}
    		
    		// If the last read doesn't have the same name, cache it
    		if (!read.getReadName().equals(readName)) {
    			cachedRead = read;
    		}
    	}
    	
    	return reads;
    }
    
    @Override
	public Iterator<List<SAMRecord>> iterator() {
		return new SamMultiMappingIterator(this);
	}
    
    public void close() {
    	inputSam.close();
    }

	private static class SamMultiMappingIterator implements Iterator<List<SAMRecord>> {

        private SamMultiMappingReader reader;
        
        SamMultiMappingIterator(SamMultiMappingReader reader) {
            this.reader = reader;
        }
        
        @Override
        public boolean hasNext() {
        	return reader.hasMoreReads();        	
        }
        
        @Override
        public List<SAMRecord> next() {
        	return reader.getNextReadList();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException("Remove not supported for SamMultiMappingIterator.");
        }
    }

}
