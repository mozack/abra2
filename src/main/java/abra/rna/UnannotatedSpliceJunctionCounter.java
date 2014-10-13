package abra.rna;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import abra.ReadBlock;

import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Counts splice junctions in a bam or sam file
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class UnannotatedSpliceJunctionCounter {
    
    private Map<SpliceJunction, Integer> spliceJunctionCounts = new HashMap<SpliceJunction, Integer>(); 
    
    
    public void count(String inputFile, String outputFile) throws IOException {
        
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFile, false));
        
        File file = new File(inputFile);
        
        SAMFileReader inputSam = new SAMFileReader(file);
        inputSam.setValidationStringency(ValidationStringency.SILENT);
        
        int count = 0;

        for (SAMRecord read : inputSam) {
            updateJunctionCounts(read);
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " reads.");
            }
        }
         
        inputSam.close();
        outputCounts(writer);
        writer.close();
    }
    
    private void outputCounts(BufferedWriter writer) throws IOException {
        System.out.println("Writing counts.");
        
        List<SpliceJunction> junctions = new ArrayList<SpliceJunction>(spliceJunctionCounts.keySet());
        
        
        for (SpliceJunction junction : junctions) {
            
            String line = junction + "\n";

            writer.write(line);
        }
    }
    
    private void updateCount(SpliceJunction junction) {
        Integer count = spliceJunctionCounts.get(junction);
        if (count == null) {
            spliceJunctionCounts.put(junction, 1);
        } else {
            spliceJunctionCounts.put(junction, count+1);
        }
    }
    
    private void updateJunctionCounts(SAMRecord read) {
        
        if (read.getCigarString().contains("N")) {
            for (ReadBlock block : ReadBlock.getReadBlocks(read)) {
                if (block.getType() == CigarOperator.N) {
                    updateCount(new SpliceJunction(read.getReferenceName(), block.getReferenceStart(), block.getReferenceStop()));
                }
            }
        }
    }
    
    public static void run(String[] args) throws IOException {
    	String input = args[0];
    	String output = args[1];
    	UnannotatedSpliceJunctionCounter counter = new UnannotatedSpliceJunctionCounter();
    	counter.count(input, output);
    }
    
    public static void main(String[] args) throws Exception {
    	run(args);
    	
    	
//        SpliceJunctionMap map = new SpliceJunctionMap("/home/lisle/gaf/splice_junctions.txt");
//        UnannotatedSpliceJunctionCounter counter = new UnannotatedSpliceJunctionCounter();
        
//        counter.count("/home/lisle/data/junction/small_sorted_by_read.bam", "/home/lisle/data/junction/counts_new.txt");
        
//        counter.count("/home/lisle/data/junction/1.sam", "/home/lisle/data/junction/counts.txt");
        
        
//        SpliceJunctionMap map = new SpliceJunctionMap(args[0]);
//        SpliceJunctionCounter counter = new SpliceJunctionCounter(map);
//        
//        counter.count(args[1], args[2]); 
    }
    
    static class SpliceJunction implements Comparable<SpliceJunction> {
    	
    	String chromosome;
    	int start;
    	int stop;
    	
    	SpliceJunction(String chr, int start, int stop) {
    		this.chromosome = chr;
    		this.start = start;
    		this.stop = stop;
    	}
    	
        @Override
        public boolean equals(Object obj) {
            SpliceJunction that = (SpliceJunction) obj;
            
            return this.chromosome.equals(that.chromosome) && this.start == that.start && this.stop == that.stop;
        }
        
        @Override
        public int hashCode() {
            int result = 17;
            
            result = 31 * result + chromosome.hashCode();
            result = 31 * result + start;
            result = 31 * result + stop;
            
            return result;
        }

        @Override
        public int compareTo(SpliceJunction that) {
            int compare = this.chromosome.compareTo(that.chromosome);
            
            if (compare == 0) {
                compare = this.start - that.start;
            }
            
            if (compare == 0) {
            	compare = this.stop - that.stop;
            }
            
            return compare;
        }
        
        public String toString() {
            return chromosome + "\t" + start + "\t" + stop;
        }
    }
}
