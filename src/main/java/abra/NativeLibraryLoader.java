/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
package abra;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;

/**
 * Manages loading of native assembly code.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class NativeLibraryLoader {
	
	public static final String ABRA     = "libAbra.so";
	public static final String SSW      = "libssw.so";
	public static final String SSW_JNI  = "libsswjni.so";
	public static final String DEFLATOR = "libIntelDeflater.so";
	
	public void load(String tempDir, String library, boolean isLenient) {
		try {
			String urlPath = "/" + library;
			
			URL url = NativeLibraryLoader.class.getResource(urlPath);
			
			if (url != null) {
				File file = new File(tempDir + "/" + library);
				
		        final InputStream in = url.openStream();
		        final OutputStream out = new BufferedOutputStream(new FileOutputStream(file));
		
		        int len = 0;
		        byte[] buffer = new byte[8192];
		        while ((len = in.read(buffer)) > -1) {
		            out.write(buffer, 0, len);
		        }
		        
		        out.close();
		        in.close();
				
				Logger.info("Loading native library from: " + file.getAbsolutePath());
				System.load(file.getAbsolutePath());
				
				file.deleteOnExit();
				
			} else if (!isLenient) {
				throw new RuntimeException("Unable to load library: " + library + " from path [" + urlPath + "] into tempdir: [" + tempDir + "]");
			}
		} catch (Throwable t) {
			Logger.error("Error loading: " + library + " from : " + tempDir);
			t.printStackTrace();
			
			if (!isLenient) {
				throw new RuntimeException(t);
			}
		}
	}
}
