package abra;

import static abra.Logger.log;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Manages threading
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ThreadManager {
	private int numThreads;
	private List<AbraRunnable> threads = new ArrayList<AbraRunnable>();

	public ThreadManager(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public Thread spawnThread(AbraRunnable runnable) {
		try {
			waitForAvailableThread();
		} catch (InterruptedException e) {}
		
		runnable.setSpawnStartTime(System.nanoTime());
		addThread(runnable);
		Thread thread = new Thread(runnable);
		thread.start();
		
		return thread;
	}
	
	private synchronized void addThread(AbraRunnable thread) {
		threads.add(thread);
	}
	
	public synchronized void removeThread(AbraRunnable thread) {
		threads.remove(thread);
	}
	
	private synchronized int activeThreads() {
		return threads.size();
	}
	
	private void waitForAvailableThread() throws InterruptedException {
		while (activeThreads() == numThreads) {
			Thread.sleep(50);
		}
	}
	
	public void waitForAllThreadsToComplete() throws InterruptedException, IOException {
		long start = System.currentTimeMillis();
		while (activeThreads() > 0) {
			long elapsedSecs = (System.currentTimeMillis() - start) / 1000;
			if ((elapsedSecs % 60) == 0) {
				log("Waiting on " + threads.size() + " threads.");
			}
			Thread.sleep(500);
		}
	}
}
