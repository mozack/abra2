/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang.ArrayUtils;

/**
 * Utility class for reversing and complementing bases.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReverseComplementor {
    
    private static final Map<Byte, Byte> complementMap = new HashMap<Byte, Byte>();
    
    static {
        complementMap.put((byte) 'C', (byte) 'G');
        complementMap.put((byte) 'G', (byte) 'C');
        complementMap.put((byte) 'T', (byte) 'A');
        complementMap.put((byte) 'A', (byte) 'T');
    }

    /**
     * Returns a new byte array containing the contents of the input byte array
     * reversed.  The input byte array is not modified.
     */
    public byte[] reverse(byte[] input) {
        byte[] bytes = ArrayUtils.clone(input);
        ArrayUtils.reverse(bytes);
        
        return bytes;
    }
   
    /**
     * Modifies the input byte array of bases to contain the complement.
     */
    public void complementInPlace(byte[] bytes) {
        for (int i=0; i<bytes.length; i++) {
            if (complementMap.containsKey(bytes[i])) {
                bytes[i] = complementMap.get(bytes[i]);
            }
        }
    }
    
    /**
     * Returns a new byte array containing the reverse complement bases of the
     * input byte array.  The input byte array is not modified.
     */
    public byte[] reverseComplement(byte[] input) {
        byte[] bytes = reverse(input);
        complementInPlace(bytes);
        
        return bytes;
    }

    /**
     * Returns the reverse complement of the input string.
     */
    public String reverseComplement(String input) {
    	return new String(reverseComplement(input.getBytes()));
    }
    
    /**
     * Returns the reverse of the input string.
     */
    public String reverse(String input) {
    	return new String(reverse(input.getBytes()));
    }
    
    public static void main(String[] args) {
    	ReverseComplementor rc = new ReverseComplementor();
    	
//    	String rcStr = rc.reverseComplement("AATGATACGGCGACCACC");
//    	String rcStr = rc.reverseComplement("AATGATACGGCGACCACCGAGATCTAC");
    	String rcStr = rc.reverseComplement("GGCTCCCAGACAGGACAGCCATCATCACAAAAGGAGATATCAAGAGGA");
    	
    	
//    	AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
    	
//		TAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT
    	
    	
    	System.err.println(rcStr);
    	
//    	rcStr = rc.reverseComplement("GACCCCCGTGACAGCTCGGGCTCC");
//    	
//    	System.out.println(rcStr);

    	
//    	String s = "ccagctatgcagagctataagtcaattaaatctgttttatttataaattac";
//    	System.out.println(s.toUpperCase());
    }
}
