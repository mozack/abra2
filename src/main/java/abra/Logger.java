package abra;

import java.util.Date;

/**
 * Simple logging class.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Logger {
	public static void log(String message) {
		System.out.println(new Date() + " : " + message);
	}
}
