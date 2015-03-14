package abra;

import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class ConcurrentQueue<T> extends ConcurrentLinkedQueue<T> {
	private AtomicInteger size = new AtomicInteger(0);
	
	@Override
	public int size() {
		return size.get();
	}
	
	public T poll() {
		size.decrementAndGet();
		return super.poll();
	}
	
	public boolean add(T t) {
		size.incrementAndGet();
		return super.add(t);
	}
}
