package abra;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class Draino implements Runnable {

	private InputStream is;
	private OutputStream os;
	
	public Draino(InputStream is, OutputStream os) {
		this.is = is;
		this.os = os;
	}
	
	@Override
	public void run() {
		try {

			byte[] buffer = new byte[4096];
			
			int bytesRead = is.read(buffer);
			
			while (bytesRead >= 0) {
				os.write(buffer, 0, bytesRead);
				bytesRead = is.read(buffer);
			}
			
			os.close();
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

}
