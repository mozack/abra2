package abra;

public class Pair<T, Y> {
	private T t;
	private Y y;
	public Pair(T t, Y y) {
		this.t = t;
		this.y = y;
	}
	
	public T getFirst() {
		return t;
	}
	
	public Y getSecond() {
		return y;
	}

}
