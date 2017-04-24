package abra;

import java.util.Date;
import java.util.HashMap;
import java.util.Map;

/**
 * Simple logging class.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class Logger {
	
	enum Level { TRACE, DEBUG, INFO, WARN, ERROR };
	
	public static Level LEVEL = Level.INFO;
	
	private static Map<String, Level> stringToLevel;
	
	static {
		stringToLevel = new HashMap<String, Level>();
		stringToLevel.put("TRACE", Level.TRACE);
		stringToLevel.put("DEBUG", Level.DEBUG);
		stringToLevel.put("INFO", Level.INFO);
		stringToLevel.put("WARN", Level.WARN);
		stringToLevel.put("ERROR", Level.ERROR);
	}
	
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
	
	public static void info(String format, Object... args) {
		if (LEVEL == Level.TRACE || LEVEL == Level.DEBUG || LEVEL == Level.INFO) {
			log(String.format(format, args), Level.INFO);
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
	
	public static void setLevel(String str) {
		Level level = stringToLevel.get(str.toUpperCase());
		if (level == null) {
			throw new IllegalArgumentException("Log level must be one of trace, debug, info, warn or error.");
		}
		
		Logger.LEVEL = level;
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
