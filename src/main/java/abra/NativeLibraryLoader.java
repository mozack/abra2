/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;

/**
 * Manages loading of native assembly code.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class NativeLibraryLoader {
	
	public void load(String tempDir, String library, boolean isLenient) {
		//String urlPath = "/libAbra.so";
		String urlPath = "/" + library;
		
		URL url = NativeLibraryLoader.class.getResource(urlPath);
		
		if (url != null) {
//			File file = new File(tempDir + "/libAbra.so");
			File file = new File(tempDir + "/" + library);
			
			try {
		        final InputStream in = url.openStream();
		        final OutputStream out = new BufferedOutputStream(new FileOutputStream(file));
		
		        int len = 0;
		        byte[] buffer = new byte[8192];
		        while ((len = in.read(buffer)) > -1) {
		            out.write(buffer, 0, len);
		        }
		        
		        out.close();
		        in.close();
			} catch (IOException e) {
				e.printStackTrace();
				throw new RuntimeException(e);
			}
			
			System.err.println("Loading native library from: " + file.getAbsolutePath());
			System.load(file.getAbsolutePath());
		} else if (!isLenient) {
			throw new RuntimeException("Unable to load library: " + library + " from path [" + urlPath + "] into tempdir: [" + tempDir + "]");
		}
	}
	
}
