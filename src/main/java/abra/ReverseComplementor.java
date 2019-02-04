/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import org.apache.commons.lang.ArrayUtils;

/**
 * Utility class for reversing and complementing bases.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ReverseComplementor {

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
     * Returns the reverse complement of the input string, non-DNA characters are allowed and just reversed.
     */
    public static String reverseComplement(String s) {
        char[] reverse = new char[s.length()];
        for (int i = 0; i < reverse.length; i++) {
            switch (s.charAt(i)) {
                case 'A': reverse[reverse.length-i-1] = 'T';break;
                case 'T': reverse[reverse.length-i-1] = 'A';break;
                case 'C': reverse[reverse.length-i-1] = 'G';break;
                case 'G': reverse[reverse.length-i-1] = 'C';break;
                default: //non-DNA input, just reverse char that was there
                    reverse[reverse.length-i-1] = s.charAt(i);
            }
        }
        return new String(reverse);
    }
    
    /**
     * Returns the reverse of the input string.
     */
    public String reverse(String input) {
    	return new String(reverse(input.getBytes()));
    }
}
