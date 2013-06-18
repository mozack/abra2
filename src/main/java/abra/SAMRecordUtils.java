package abra;

import java.util.ArrayList;
import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

public class SAMRecordUtils {

	
	public static void removeHardClips(SAMRecord read) {
		Cigar cigar = read.getCigar();
		
		CigarElement firstElement = cigar.getCigarElement(0);
		CigarElement lastElement  = cigar.getCigarElement(cigar.numCigarElements()-1);
		
		if ((firstElement.getOperator() == CigarOperator.H) ||
			(lastElement.getOperator() == CigarOperator.H)) {
			
			Cigar newCigar = new Cigar();
			for (CigarElement element : cigar.getCigarElements()) {
				if (element.getOperator() != CigarOperator.H) {
					newCigar.add(element);
				}
			}
			
			read.setCigar(newCigar);
		}
	}
	
	/**
	 * Remove leading or trailing soft clips from the input read.
	 * Does not modify a read entirely comprised of soft clips.
	 */
	public static void removeSoftClips(SAMRecord read) {
		
		Cigar cigar = read.getCigar();
		
		CigarElement firstElement = cigar.getCigarElement(0);
		CigarElement lastElement  = cigar.getCigarElement(cigar.numCigarElements()-1);
		
		if ((firstElement.getOperator() == CigarOperator.S) ||
			(lastElement.getOperator() == CigarOperator.S)) {
		
			Cigar newCigar = new Cigar();
			
			String bases = read.getReadString();
			//String qualities = read.getBaseQualityString();
					
			if (firstElement.getOperator() == CigarOperator.S) {
				bases = bases.substring(firstElement.getLength(), bases.length());
				//qualities = qualities.substring(firstElement.getLength(), qualities.length()-1);
			} else {
				newCigar.add(firstElement);
			}
			
			for (int i=1; i<cigar.numCigarElements()-1; i++) {
				newCigar.add(cigar.getCigarElement(i));
			}
			
			if (lastElement.getOperator() == CigarOperator.S) {
				bases = bases.substring(0, bases.length() - lastElement.getLength());
				//qualities = qualities.substring(0, qualities.length() - lastElement.getLength() - 1);
			} else {
				newCigar.add(lastElement);
			}
			
			read.setCigar(newCigar);
			read.setReadString(bases);
			//read.setBaseQualityString(qualities);
		}
	}

	public static Cigar subset(Cigar cigar, int startIdx, int endIdx) {
		
		List<CigarElement> subElements = new ArrayList<CigarElement>();
		
		// Find first element and index into first element
		int len = 0;
		int elemIdx = -1;
		for (CigarElement elem : cigar.getCigarElements()) {

			// Treat deletions as zero length.
			int elemLength = elem.getOperator() == CigarOperator.D ? 0 : elem.getLength();
			
			if (elemIdx < 0) {
				
				// Find first element (Should never be a deletion)
				int elemStart = len;
				int elemStop  = len + elemLength;
				
				if ((startIdx >= elemStart) && (startIdx < elemStop)) {
					elemIdx = startIdx - elemStart;
					int elemLen = Math.min(elem.getLength()-elemIdx, endIdx-startIdx+1);
					CigarElement newElem = new CigarElement(elemLen, elem.getOperator());
					subElements.add(newElem);
				}
				
				len += elemLength;
				
			} else if ((len + elemLength) <= endIdx) {
				// Add this entire element
				subElements.add(elem);
				len += elemLength;
			} else if (len <= endIdx) {
				// Add part of last sub element (should never be a deletion)
				CigarElement newElem = new CigarElement(endIdx-len+1, elem.getOperator());
				subElements.add(newElem);
				break;
			} else {
				break;
			}
		}
		
		return new Cigar(subElements);
	}
}
