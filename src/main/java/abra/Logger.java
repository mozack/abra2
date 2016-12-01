package abra;

import java.util.Date;

/**
 * Simple logging class.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Logger {
	
	enum Level { TRACE, DEBUG, INFO, WARN, ERROR };
	
	public static Level LEVEL = Level.INFO;
	
	// For trace or debug messages, use varargs to avoid string concatenation unless enabled
	
	public static void trace(String format, Object... args) {
		if (LEVEL == Level.TRACE) {
			log(String.format(format, args), Level.TRACE);
		}
	}
	
	public static void debug(String format, Object... args) {
		if (LEVEL == Level.DEBUG || LEVEL == Level.TRACE) {
			log(String.format(format, args), Level.DEBUG);
		}
	}
	
	public static void info(String message) {
		if (LEVEL == Level.TRACE || LEVEL == Level.DEBUG || LEVEL == Level.INFO) {
			log(message, Level.INFO);
		}
	}
	
	public static void warn(String message) {
		if (LEVEL == Level.TRACE || LEVEL == Level.DEBUG || LEVEL == Level.INFO || LEVEL == Level.WARN) {
			log(message, Level.WARN);
		}
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
			case TRACE:
				levelStr = "TRACE";
				break;
		}
		
		System.err.println(levelStr + "\t" + new Date() + "\t" + message);
	}
}
