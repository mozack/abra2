package abra;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Utility class for outputting fastq files
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class FastqOutputFile {
    private BufferedWriter writer;
    
    public void init(String filename) throws IOException {
        writer = new BufferedWriter(new FileWriter(filename, false));
    }
    
    public void write(FastqRecord record) throws IOException {
        for (int i = 0; i<record.getLines().length; i++) {
            writer.write(record.getLines()[i]);
            writer.write('\n');
        }
    }
    
    public void close() throws IOException {
        writer.close();
    }
}

