package abra.cadabra;

import org.apache.commons.math3.special.Gamma;

public class BetaBinomial {
	
	private static double fita = 0.4657371;
	private static double fitb = 1494.0276936;
	private static double minVal = 1e-16;
		
	public static double betabinCDF(int depth, double maxVal) {
        double v = 0;
        for (int i = depth; i >= maxVal; i--) {
            double x = betabinPMFG(depth, i, fita, fitb);
            v += x;
        }
        return pinch(v);
    }
	
    private static double pinch(double n) {
        double r = n;
        if (r <= minVal) {
            return minVal;
        }
        if (r >= 1 - minVal) {
            r = 1 - minVal;
            return r;
        }
        return r;
    }
	
    private static double betabinPMFG(int n, int k, double a, double b) {
        double b1 = Gamma.logGamma(n + 1);
        double b2 = Gamma.logGamma(k + 1) + Gamma.logGamma(n - k + 1);
        double b3 = Gamma.logGamma(k + a) + Gamma.logGamma(n - k + b);
        double b4 = Gamma.logGamma(n + a + b);
        double b5 = Gamma.logGamma(a + b);
        double b6 = Gamma.logGamma(a) + Gamma.logGamma(b);
        double v = b1 - b2 + b3 - b4 + b5 - b6;
        return Math.exp(v);
    }
}
