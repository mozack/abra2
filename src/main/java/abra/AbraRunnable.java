package abra;

/**
 * Abstract base class for ABRA threads.
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public abstract class AbraRunnable implements Runnable {

	long spawnStartTime;
	private ThreadManager threadManager;
	
	public AbraRunnable(ThreadManager threadManager) {
		this.threadManager = threadManager;
	}
	
	@Override
	public void run() {
		try {
			go();
		} catch (Throwable t) {
			t.printStackTrace();
			System.exit(-1);
		} finally {
			threadManager.removeThread(this);
		}
	}
	
	public void setSpawnStartTime(long spawnStartTime) {
		this.spawnStartTime = spawnStartTime;
	}
	
	public abstract void go() throws Exception;
}
