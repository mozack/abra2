package abra;

public class MutableBoolean {

	private boolean value = false;
	
	public boolean isTrue() {
		return value;
	}
	
	public void setValue(boolean val) {
		this.value = val;
	}
}
