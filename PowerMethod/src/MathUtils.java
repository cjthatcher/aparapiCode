
/*
 * @Author Chris Thatcher Jan 16. 2014
 */
public class MathUtils {
	
	/* Computes the relative error between two doubles.
	 * 
	 * That is, returns |newValue - oldValue| / |oldValue|
	 * 
	 * @param double the old value
	 * @param double the new value
	 * @return the relative error between the two values
	 */
	public static float computeRelativeError(float oldValue, float newValue)
	{
		return Math.abs(newValue - oldValue) / Math.abs(oldValue);
	}
	
	
	
	
	
	
	

	
}
