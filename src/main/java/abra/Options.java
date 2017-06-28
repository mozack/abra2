/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.IOException;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * Abstract base class for helping with options parsing.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public abstract class Options {
	protected static final String HELP = "help";
	
	private OptionSet options;

    protected void printHelp() {
        try {
        	getOptionParser().printHelpOn(System.err);
        }
        catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException("IOException encountered when attempting to output help.");
        }
    }
    
    public void parseOptions(String[] args) {

        try {
            options = getOptionParser().parse(args);
        
            if (options.has(HELP)) {
                printHelp();
            } else {
            	init();
                validate();
            }
        } catch (joptsimple.OptionException e) {
            System.err.println(e.getMessage());
            printHelp();
            throw e;
        }
    }
    
    protected OptionSet getOptions() {
    	return options;
    }

    abstract protected OptionParser getOptionParser();
    
    abstract protected void validate();
    
    protected void init() {
    }
}
