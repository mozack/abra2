package abra;

import static abra.Logger.log;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Manages threading
 * 
 * @author Lisle E. Mose (lmose at unc dot edu)
 */
public class ThreadManager {
	
	private static final int MAX_PENDING = 100;
	
	private int numThreads;
	private List<AbraRunnable> threads = new ArrayList<AbraRunnable>();
	private ExecutorService executor;

	public ThreadManager(int numThreads) {
		this.numThreads = numThreads;
		executor = Executors.newFixedThreadPool(numThreads);
	}
	
	public void spawnThread(AbraRunnable runnable) {
		
		try {
			waitForAvailableThread();
		} catch (InterruptedException e) {}
		
		addThread(runnable);
		
		executor.submit(runnable);
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
		while (activeThreads() >= MAX_PENDING) {
			Thread.sleep(50);
		}
	}
	
	public void waitForAllThreadsToComplete() throws InterruptedException, IOException {
		executor.shutdown();
		while (!executor.awaitTermination(60, TimeUnit.SECONDS)) {
			Logger.info("Waiting on " + threads.size() + " threads.");
		}
	}

	public int getNumThreads() {
		return numThreads;
	}
}
