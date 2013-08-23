/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.util.Arrays;

public class Sequence {
	private short length;
	private byte[] sequence;
//	private int hashCode;
	
	private static final byte A = 0;
	private static final byte C = 1;
	private static final byte T = 2;
	private static final byte G = 3;
	
	public Sequence(String str) {
		length = (short) str.length();
		
		int numBits  = str.length() * 2;
		int numBytes = numBits / 8;
		if ((numBits % 8) > 0) {
			numBytes += 1;
		}
		
		sequence = new byte[numBytes];
		Arrays.fill(sequence, (byte) 0);
		
		for (short i=0; i<str.length(); i++) {
			char ch = str.charAt(i);
			byte base = getBase(str, ch);
			
			short bucketIdx = (short) (i * 2 / 8);
			
			byte bucket = sequence[bucketIdx];
			
			short posInBucket = (short) ((i*2) % 8);
			
			byte positionedBase = (byte) (base << posInBucket);
			
			bucket = (byte) (bucket | positionedBase);
			
			sequence[bucketIdx] = bucket;
		}
	}
	
	public char getFirstCharacter() {
		byte bucket = sequence[0];
		byte base = (byte) (bucket & 3);
		return getChar(base);
	}
	
	public String getSequenceAsString() {
		StringBuffer seq = new StringBuffer();
		
		for (int i=0; i<length; i++) {
			short bucketIdx = (short) (i * 2 / 8);
			
			byte bucket = sequence[bucketIdx];
			
			short posInBucket = (short) ((i*2) % 8);
			
			bucket = (byte) (bucket >> posInBucket);
			
			// (filter all but last 2 bits) i.e. bucket & 0x00000011
			byte base = (byte) (bucket & 3);
			
			char ch = getChar(base);
			
			seq.append(ch);
		}
		
		return seq.toString();
	}
	
	public int hashCode() {
		return Arrays.hashCode(sequence);
	}
	
	public boolean equals(Object obj) {
		Sequence that = (Sequence) obj;
		return Arrays.equals(this.sequence, that.sequence);
	}
	
	private char getChar(byte base) {
		switch (base) {
			case A:
				return 'A';
			case C:
				return 'C';
			case T:
				return 'T';			
			case G:
				return 'G';
			default:
				throw new IllegalArgumentException("Invalid base: " + base);
		}
	}
	
	private byte getBase(String str, char ch) {
		switch (ch) {
			case 'A':
				return A;
			case 'C':
				return C;
			case 'T':
				return T;			
			case 'G':
				return G;
			default:
				throw new IllegalArgumentException("Invalid base: " + ch + " for sequence: " + str);
		}
	}
	
	public static void main(String[] args) {
		
		/*
		int i = 2;
		
		i = i >> 1;
		
		System.out.println(Integer.toBinaryString(i));
		
//		int i = 2 >> 2;
		
		System.out.println(i);
		*/
		
		Sequence s1 = new Sequence("ATCGATCG");
		Sequence s2 = new Sequence("ATCGATCC");
		Sequence s3 = new Sequence("ATCGATCG");
		Sequence s4 = new Sequence("ATCGATCGG");
		Sequence s5 = new Sequence("TCGATCGG");
		
		System.out.println(s1.getSequenceAsString());
		System.out.println(s2.getSequenceAsString());
		System.out.println(s3.getSequenceAsString());
		System.out.println(s4.getSequenceAsString());
		System.out.println(s5.getSequenceAsString());
		
		System.out.println(s1.equals(s2));
		System.out.println(s2.equals(s1));
		System.out.println(s1.equals(s3));
		System.out.println(s3.equals(s1));
		System.out.println(s3.equals(s4));
		System.out.println(s3.equals(s5));
		System.out.println(s5.equals(s5));
		
		
		System.out.println(s1.getFirstCharacter());
		System.out.println(s5.getFirstCharacter());
	}
}
