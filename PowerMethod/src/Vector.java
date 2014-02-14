/*
 * The Vector class is a subclass of Matrix. It is a Matrix that only has one row or column.
 * 
 * It's row/column is held in a variable called "data" and can get fetched using the getData() call.
 * 
 * 
 */

/*
 * @Author Chris Thatcher Jan 16. 2014
 */
public class Vector extends Matrix {
	private final Orientation orientation;
	private float[] data;
	
	/*
	 * Static factory method that returns a Vector given one array of floats and its orientation.
	 * Constructs a Vector from one array of floats. The Orientation enum is needed
	 * to decide whether this is a row or column vector.
	 * 
	 * @param float[] An array holding all the elements of the vector
	 * @param Orientation The orientation of the vector (whether it is row or column)
	 * 
	 * @return A new vector with the correct orientation and data
	 */
	public static Vector makeVector(float[] column, Orientation orientation)
	{
		return new Vector(column, orientation);
	}
	
	/*
	 * Constructs a Vector from one array of floats. The Orientation enum is needed
	 * to decide whether this is a row or column vector.
	 * 
	 * @param float[] An array holding all the elements of the vector
	 * @param Orientation The orientation of the vector (whether it is row or column)
	 */
	public Vector(float[] column, Orientation orientation)
	{
		super(column, orientation);
		this.orientation = orientation;
		this.data = column;
	}
	
	/*
	 * Constructs a vector from a 2D array of values. The orientation of the vector is 
	 * inferred based on the length of the arrays passed in. Note that this constructor is unsafe
	 * and should really not be used.
	 * 
	 * @param float[][] the underlying matrix data. This should either be one row or one column
	 */
	public Vector(float[][] values) {
		super(values);
		if (values.length > values[0].length)
		{
			this.orientation = Orientation.ROW;
			this.data = values[0];
		}
		else
		{
			this.orientation = Orientation.COLUMN;
			
			this.data = new float[values.length];
			for (int i = 0; i < values.length; i++)
			{
				data[i] = values[0][i];
			}
		}
	}
	

	/*
	 * Returns the transpose of the current Vector. Not that the original vector is unchanged. This method returns a free-standing copy of the original
	 * vector with the opposite orientation. 
	 * 
	 * @return A new Vector that is equivalent to the calling Vector's transpose.
	 */
	public Vector transpose() {
		//We need to make a clone of the original data.
		float[] f = data.clone();
		return new Vector(f, orientation == Orientation.ROW ? Orientation.COLUMN : Orientation.ROW);
	}
	
	/*
	 * Getter for orientation.
	 * 
	 * @return The orientation of the calling Vector
	 */
	public Orientation getOrientation(){
		return orientation;
	}
	
	/*
	 * Getter for data. This returns an array of floats that represents the row (or column) of the vector.
	 * 
	 * @return The row (or column) of elements in the Vector as an array of floats.
	 */
	public float[] getData() {
		return data;
	}
	
	/*
	 * Computes the magnitude of a given calling Vector v. The magnitude is found
	 * by taking the square root of the sum of squares of each element in the vector.
	 * 
	 * @return the magnitude of the calling vector as a float
	 */
	public float getMagnitude()
	{
		float[] values = this.getData();
		float sumOfSquares = 0;
		
		for (float d : values)
		{
			sumOfSquares += Math.pow(d, 2);
		}
		
		return (float) Math.pow(sumOfSquares, 0.5);
	}
	
	/*
	 * Returns the normalized form of the calling vector. The magnitude of the 
	 * vector is found, and then a new vector scaled by 1 / magnitude is returned.
	 * The original vector remains unchanged. Note that the magnitude of the 
	 * resultant vector is 1.
	 * 
	 * @return a unit vector in the direction of the calling vector
	 */
	public Vector getNormalizedVector()
	{
		float magnitude = getMagnitude();
		
		float[] oldData = this.getData();
		float[] nuevoData = new float[oldData.length];

		for (int i = 0; i < oldData.length; i++)
		{
			nuevoData[i] = oldData[i] / magnitude;
		}
		
		return new Vector(nuevoData, this.getOrientation());
	}
}
