package abra;

import java.util.Date;

/**
 * Simple logging class.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Logger {
	
	enum Level { DEBUG, INFO, WARN, ERROR };
	
	public static Level LEVEL = Level.INFO;
	
	// For debug messages, use varargs to avoid string concatenation unless necessary
	public static void debug(Object ... parts) {
		if (LEVEL == Level.DEBUG) {
			StringBuffer msg = new StringBuffer();
			for (Object part : parts) {
				msg.append(part);
			}
			log(msg.toString(), Level.DEBUG);
		}
	}
	
	public static void info(String message) {
		log(message, Level.INFO);
	}
	
	public static void warn(String message) {
		log(message, Level.WARN);
	}
	
	public static void error(String message) {
		log(message, Level.ERROR);
	}
	
	public static void log(String message, Level level) {
		
		String levelStr = "UNKNOWN";
		
		switch (level) {
			case ERROR:
				levelStr = "ERROR";
				break;
			case WARN:
				levelStr = "WARNING";
				break;
			case INFO:
				levelStr = "INFO";
				break;
			case DEBUG:
				levelStr = "DEBUG";
				break;
		}
		
		System.err.println(levelStr + "\t" + new Date() + "\t" + message);
	}
}
