/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;


import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;

import static org.testng.Assert.assertEquals;

import org.testng.Assert;
import org.testng.annotations.Test;

public class IndelShifterTest {
	
	private IndelShifter indelShifter = new IndelShifter();

	@Test (groups = "unit" )
	public void testShiftCigarLeft_basic() {
		Cigar cigar = new Cigar();
		
		cigar.add(new CigarElement(10, CigarOperator.M));
		cigar.add(new CigarElement(3, CigarOperator.D));
		cigar.add(new CigarElement(40, CigarOperator.M));
		
		Cigar newCigar;

		newCigar = indelShifter.shiftCigarLeft(cigar, 10);
		Assert.assertEquals(newCigar.toString(), "3D50M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 9);
		Assert.assertEquals(newCigar.toString(), "1M3D49M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 8);
		Assert.assertEquals(newCigar.toString(), "2M3D48M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 4);
		Assert.assertEquals(newCigar.toString(), "6M3D44M");

		newCigar = indelShifter.shiftCigarLeft(cigar, 2);
		Assert.assertEquals(newCigar.toString(), "8M3D42M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 1);
		Assert.assertEquals(newCigar.toString(), "9M3D41M");
	}
	
	@Test (groups = "unit" )
	public void testShiftCigarLeft_softClipping() {
		Cigar cigar = new Cigar();
		
		cigar.add(new CigarElement(2, CigarOperator.S));
		cigar.add(new CigarElement(6, CigarOperator.M));
		cigar.add(new CigarElement(2, CigarOperator.I));
		cigar.add(new CigarElement(30, CigarOperator.M));
		cigar.add(new CigarElement(10, CigarOperator.S));
		
		Cigar newCigar;

		newCigar = indelShifter.shiftCigarLeft(cigar, 6);
		Assert.assertEquals(newCigar.toString(), "2S2I36M10S");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 5);
		Assert.assertEquals(newCigar.toString(), "2S1M2I35M10S");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 4);
		Assert.assertEquals(newCigar.toString(), "2S2M2I34M10S");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 3);
		Assert.assertEquals(newCigar.toString(), "2S3M2I33M10S");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 2);
		Assert.assertEquals(newCigar.toString(), "2S4M2I32M10S");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 1);
		Assert.assertEquals(newCigar.toString(), "2S5M2I31M10S");
	}
	
	@Test (groups = "unit" )
	public void testShiftCigarLeft_insertAtTail() {
		Cigar cigar = new Cigar();
		
		cigar.add(new CigarElement(40, CigarOperator.M));
		cigar.add(new CigarElement(10, CigarOperator.I));
		
		Cigar newCigar;

		newCigar = indelShifter.shiftCigarLeft(cigar, 40);
		Assert.assertEquals(newCigar.toString(), "10I40M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 39);
		Assert.assertEquals(newCigar.toString(), "1M10I39M");

		newCigar = indelShifter.shiftCigarLeft(cigar, 30);
		Assert.assertEquals(newCigar.toString(), "10M10I30M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 1);
		Assert.assertEquals(newCigar.toString(), "39M10I1M");
	}
	
	@Test (groups = "unit" )
	public void testShiftCigarLeft_multipleIndels() {
		Cigar cigar = new Cigar();
		
		cigar.add(new CigarElement(20, CigarOperator.M));
		cigar.add(new CigarElement(1, CigarOperator.I));
		cigar.add(new CigarElement(5, CigarOperator.M));
		cigar.add(new CigarElement(3, CigarOperator.D));
		cigar.add(new CigarElement(24, CigarOperator.M));
		
		Cigar newCigar;
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 20);
		Assert.assertEquals(newCigar.toString(), "1I5M3D44M");
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 10);
		Assert.assertEquals(newCigar.toString(), "10M1I5M3D34M");

		newCigar = indelShifter.shiftCigarLeft(cigar, 1);
		Assert.assertEquals(newCigar.toString(), "19M1I5M3D25M");
	}
	
	@Test (groups = "unit" )
	public void testShiftCigarLeft_complex() {
		//3S69M1I18M1D9M
		Cigar cigar = new Cigar();
		
		cigar.add(new CigarElement(3, CigarOperator.S));
		cigar.add(new CigarElement(69, CigarOperator.M));
		cigar.add(new CigarElement(1, CigarOperator.I));
		cigar.add(new CigarElement(18, CigarOperator.M));
		cigar.add(new CigarElement(1, CigarOperator.D));
		cigar.add(new CigarElement(9, CigarOperator.M));
		
		Cigar newCigar;
		
		newCigar = indelShifter.shiftCigarLeft(cigar, 1);
		assertEquals(newCigar.toString(), "3S68M1I18M1D10M");
	}
	
	@Test (groups = "unit" )
	public void testShiftDelLeft() throws Exception {
		Cigar cigar = TextCigarCodec.decode("6M2D8M");
		Cigar newCigar = indelShifter.shiftCigarLeft(cigar, 4);
		assertEquals(TextCigarCodec.encode(newCigar), "2M2D12M");
	}

	@Test (groups = "unit" )
	public void testShiftIndelsLeft() throws Exception {
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init("test-data/test.fa");
		/*
		TCGAATCGATATATTTCCGGAACAGACTCAG
		------CGATAT--TTCCGGAA--------- <-- orig
		------CG--ATATTTCCGGAA--------- <-- new
		1234567890123456789012
		*/
		
		int refStart = 7;
		int refEnd = 22;
		Cigar cigar = TextCigarCodec.decode("6M2D8M");
		String seq = "CGATATTTCCGGAA";
		
		// 1 based input
		Cigar newCigar = indelShifter.shiftIndelsLeft(refStart, refEnd, "seq1", cigar, seq, c2r);
		assertEquals(TextCigarCodec.encode(newCigar), "2M2D12M");
	}
	
	@Test (groups = "unit" )
	public void testShiftIndelsLeft_LocalRef() throws Exception {
		
		CompareToReference2 c2r = new CompareToReference2();
		c2r.initLocal("seq1", "TCGAATCGATATATTTCCGGAACAGACTCAG");
		//c2r.init("test-data/test.fa");
		/*
		TCGAATCGATATATTTCCGGAACAGACTCAG
		------CGATAT--TTCCGGAA--------- <-- orig
		------CG--ATATTTCCGGAA--------- <-- new
		1234567890123456789012
		*/
		
		int refStart = 7;
		int refEnd = 22;
		Cigar cigar = TextCigarCodec.decode("6M2D8M");
		String seq = "CGATATTTCCGGAA";
		
		// 1 based input
		Cigar newCigar = indelShifter.shiftIndelsLeft(refStart, refEnd, "seq1", cigar, seq, c2r);
		assertEquals(TextCigarCodec.encode(newCigar), "2M2D12M");
	}
}
