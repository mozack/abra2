package abra;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

/**
 * Represents a Block of a read.  i.e. A Cigar of 15M5I30M would be represented
 * by 3 ReadBlocks.
 * 
 * @author Lisle Mose (lmose at unc dot edu)
 */
public class ReadBlock {
    private int readStart;
    private int referenceStart;
    private int length;
    private CigarOperator type;
    
    ReadBlock(int readStart, int referenceStart, int length, CigarOperator type) {
        this.readStart = readStart;
        this.referenceStart = referenceStart;
        this.length = length;
        this.type = type;
    }

    public int getReadStart() {
        return readStart;
    }

    public int getReferenceStart() {
        return referenceStart;
    }
    
    public int getReferenceStop() {
    	//TODO: This doesn't appear to be correct for delete block type
        return referenceStart + length - 1;
    }

    public int getLength() {
        return length;
    }
    
    public void setLength(int length) {
    	this.length = length;
    }
    
    public void setReferenceStart(int referenceStart) {
    	this.referenceStart = referenceStart;
    }
    
    public int getReferenceLength() {
    	int refLength = 0;
    	
    	if (type != CigarOperator.D) {
    		refLength = length;
    	}
    	
    	return refLength;
    }

    public CigarOperator getType() {
        return type;
    }
    
    /**
     * Returns a ReadBlock as a subset of this readblock
     */
    public ReadBlock getSubBlock(int accumulatedLength, int positionInRead, int maxLength) {
    	
    	// zero base + zero base - 1 base  + 1 = zero base
    	int positionInBlock = positionInRead + accumulatedLength - readStart + 1;
    	
    	if ((type == CigarOperator.N) || (type == CigarOperator.D)) {
    		// Intron / Deletion: return entire block
        	return new ReadBlock(accumulatedLength+1, referenceStart + positionInBlock, length-positionInBlock, type);
    	} else if (type == CigarOperator.S) {
    		// Soft clipped blocks begin at next block's referenceStart.  No need to change it.
    		return new ReadBlock(accumulatedLength+1, referenceStart, Math.min(maxLength, length - positionInBlock), type);
    	} else {
    		return new ReadBlock(accumulatedLength+1, referenceStart + positionInBlock, Math.min(maxLength, length - positionInBlock), type);
    	}
    }
    
    public static String toCigarString(List<ReadBlock> blocks) {
    	StringBuffer cigar = new StringBuffer();
    	
    	for (ReadBlock block : blocks) {
    		cigar.append(block.getLength());
    		cigar.append(block.getType());
    	}
    	
    	return cigar.toString();
    }
    
    public static int getTotalLength(Collection<ReadBlock> blocks) {
    	int length = 0;
    	
    	for (ReadBlock block : blocks) {
    		if (block.getType() != CigarOperator.D) { 
    			length += block.getLength();
    		}
    	}
    	
    	return length;
    }
    
    //TODO - Move elsewhere and make non-static
    public static List<ReadBlock> getReadBlocks(SAMRecord read) {    
        final Cigar cigar = read.getCigar();
        if (cigar == null) return Collections.emptyList();

        final List<ReadBlock> readBlocks = new ArrayList<ReadBlock>();
        int readBase = 1;
        int refBase  = read.getAlignmentStart();

        for (final CigarElement e : cigar.getCigarElements()) {
            
            readBlocks.add(new ReadBlock(readBase, refBase, e.getLength(), e.getOperator()));
            
            int[] basePositions = updateBasePositions(readBase, refBase, e.getOperator(), e.getLength());
            readBase = basePositions[0];
            refBase = basePositions[1];
        }
        
        return Collections.unmodifiableList(readBlocks);
    }
    
    private static int[] updateBasePositions(int readBase, int refBase, CigarOperator operator, int elemLength) {
    	switch (operator) {
	//      case H : break; // ignore hard clips
	//      case P : break; // ignore pads
	      case S : readBase += elemLength;
	//      System.out.println(read.getReadName());
	      break; // soft clip read bases
	      case N : refBase += elemLength; break;  // reference skip
	      case D : refBase += elemLength; break;
	      case I : readBase += elemLength; break;
	      case M :
	//      case EQ :
	//      case X :
	          final int length = elemLength;
	          readBase += length;
	          refBase  += length;
	          break;
	      default : throw new IllegalStateException(
	              "Case statement didn't deal with cigar op: " + operator);
    	}
    	
    	return new int[] { readBase, refBase };
    }
    
    public static void fillToLength(List<ReadBlock> blocks, int readLength) {
		int cigarLength = ReadBlock.getTotalLength(blocks);
		
		if (cigarLength < readLength) {
			ReadBlock lastBlock = blocks.get(blocks.size() - 1);
			if ((lastBlock.getType() == CigarOperator.M) ||
				(lastBlock.getType() == CigarOperator.S)) {
				// Last block is M or S, so extend it.
				lastBlock.setLength(lastBlock.getLength() + readLength - cigarLength);
			} else {
				// Last block is not M or S.  Add a new M block.
				int[] basePositions = updateBasePositions(lastBlock.getReadStart(), lastBlock.getReferenceStart(), lastBlock.getType(), lastBlock.getLength());
				int readStart = basePositions[0];
				int refStart = basePositions[1];
				blocks.add(new ReadBlock(readStart, refStart, readLength-cigarLength, CigarOperator.M));
			}
		}
    }
}
