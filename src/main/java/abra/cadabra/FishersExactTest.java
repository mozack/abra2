package abra.cadabra;

import java.util.Random;

/**
 * Basic Fisher's exact test implementation.
 * 
 * Caches factorial values in log space across calculations, thus speeding things up (a bit).
 *  
 * @author lmose
 */
public class FishersExactTest {

	private static int MAX_SIZE = 5000;
	
	// Cache of factorial values in log space
	private static double[] factorialCache = new double[MAX_SIZE+1];
	
	static {
		init();
	}
	
	private static void init() {
		for (int i=1; i<=MAX_SIZE; i++) {
			factorialCache[i] = factorialCache[i-1] + Math.log(i);
		}
	}
	
	public double oneTailedTest(int normalRef, int normalAlt, int tumorRef, int tumorAlt) {
		int row1Col1 = normalRef;
		int row1Col2 = normalAlt;
		int row2Col1 = tumorRef;
		int row2Col2 = tumorAlt;
		
		int n = row1Col1 + row1Col2 + row2Col1 + row2Col2;
		if (n > MAX_SIZE) {
			double scale = (double) MAX_SIZE / (double) n;
			row1Col1 = (int) (row1Col1 * scale);
			row1Col2 = (int) (row1Col2 * scale);
			row2Col1 = (int) (row2Col1 * scale);
			row2Col2 = (int) (row2Col2 * scale);
			
			n = row1Col1 + row1Col2 + row2Col1 + row2Col2;
		}

		int row1Sum = row1Col1 + row1Col2;
		int row2Sum = row2Col1 + row2Col2;
		int col1Sum = row1Col1 + row2Col1;
		int col2Sum = row1Col2 + row2Col2;
		
		double numerator = factorialCache[row1Sum] + factorialCache[row2Sum] + factorialCache[col1Sum] + factorialCache[col2Sum];
		
		double pObserved = getPForTable(row1Col1, row1Col2, row2Col1, row2Col2, n, numerator);
		double pValue = pObserved;
		
		while (row1Col2 > 0 && row2Col1 > 0) {
			row1Col1++;
			row1Col2--;
			row2Col1--;
			row2Col2++;
			
			double nextP = getPForTable(row1Col1, row1Col2, row2Col1, row2Col2, n, numerator);
			
			if (nextP <= pObserved) {
				pValue += nextP;
			}
		}
		
		// Cap p-value at 1 to guard against rounding errors
		return Math.min(pValue, 1.0);
	}
	
	//TODO: Extract shared code
	public double twoTailedTest(int normalRef, int normalAlt, int tumorRef, int tumorAlt) {
		int row1Col1 = normalRef;
		int row1Col2 = normalAlt;
		int row2Col1 = tumorRef;
		int row2Col2 = tumorAlt;
		
		int n = row1Col1 + row1Col2 + row2Col1 + row2Col2;
		if (n > MAX_SIZE) {
			double scale = (double) MAX_SIZE / (double) n;
			row1Col1 = (int) (row1Col1 * scale);
			row1Col2 = (int) (row1Col2 * scale);
			row2Col1 = (int) (row2Col1 * scale);
			row2Col2 = (int) (row2Col2 * scale);
			
			n = row1Col1 + row1Col2 + row2Col1 + row2Col2;
		}
		
		int row1Col1Start = row1Col1;
		int row1Col2Start = row1Col2;
		int row2Col1Start = row2Col1;
		int row2Col2Start = row2Col2;

		int row1Sum = row1Col1 + row1Col2;
		int row2Sum = row2Col1 + row2Col2;
		int col1Sum = row1Col1 + row2Col1;
		int col2Sum = row1Col2 + row2Col2;
		
		double numerator = factorialCache[row1Sum] + factorialCache[row2Sum] + factorialCache[col1Sum] + factorialCache[col2Sum];
		
		double pObserved = getPForTable(row1Col1, row1Col2, row2Col1, row2Col2, n, numerator);
		double pValue = pObserved;
		
		while (row1Col2 > 0 && row2Col1 > 0) {
			row1Col1++;
			row1Col2--;
			row2Col1--;
			row2Col2++;
			
			double nextP = getPForTable(row1Col1, row1Col2, row2Col1, row2Col2, n, numerator);
			
			if (nextP <= pObserved) {
				pValue += nextP;
			}
		}
		
		// Now the other way...
		row1Col1 = row1Col1Start;
		row1Col2 = row1Col2Start;
		row2Col1 = row2Col1Start;
		row2Col2 = row2Col2Start;
		
		while (row1Col1 > 0 && row2Col2 > 0) {
			row1Col1--;
			row1Col2++;
			row2Col1++;
			row2Col2--;
			
			double nextP = getPForTable(row1Col1, row1Col2, row2Col1, row2Col2, n, numerator);
			
			if (nextP <= pObserved) {
				pValue += nextP;
			}
		}
		
		// Cap p-value at 1 to guard against rounding errors
		return Math.min(pValue, 1.0);

	}
	
	private double getPForTable(int r1c1, int r1c2, int r2c1, int r2c2, int n, double numerator) {
		//TODO: Remove this as optimization
		if ((r1c1 + r1c2 + r2c1 + r2c2) != n) throw new IllegalArgumentException("Invalid contigency table");
		
		double denominator = factorialCache[r1c1] + factorialCache[r1c2] + factorialCache[r2c1] + factorialCache[r2c2] + factorialCache[n];  
		return Math.exp(numerator - denominator);
	}
	
	private static int nextRand(Random r) {
		return r.nextInt(10000);
	}
	
	public static void main(String[] args) {
//		int nr = 1500; int na = 110; int tr = 1400; int ta = 1100;
		
		int nr = 709; int na = 20; int tr = 711; int ta = 85;
		
		FishersExactTest t = new FishersExactTest();
		
		double p = t.oneTailedTest(nr, na, tr, ta);
		
		
//		Random r = new Random();
//		
//		long s = System.currentTimeMillis();
//		
//		for (int i=0; i<10000; i++) {
//			double p = t.oneTailedTest(nr + nextRand(r), na + nextRand(r), tr + nextRand(r), ta + nextRand(r));
//			if (i%1000 == 0) {
//				System.out.println(p);
//			}
//		}
//		
//		long e = System.currentTimeMillis();
		System.out.println("p: " + p);
		System.out.println("phred: " + (-10 * Math.log10(p)));
//		
//		System.out.println(e-s);
		
		p = t.twoTailedTest(nr, na, tr, ta);
		
		System.out.println("2 tailed: " + p);
	}
}
