package abra;

public class Clock {

	private String descriptor;
	private long startMsecs;
	private long stopMsecs;
	
	public Clock(String descriptor) {
		this.descriptor = descriptor;
	}
	
	public void start() {
		this.startMsecs = System.currentTimeMillis();
	}
	
	public long elapsedSeconds() {
		return (stopMsecs - startMsecs) / 1000;
	}
	
	public void stopAndPrint() {
		this.stopMsecs = System.currentTimeMillis();
		
		System.out.println("Clock time in " + descriptor + ": " + elapsedSeconds());
	}
}
