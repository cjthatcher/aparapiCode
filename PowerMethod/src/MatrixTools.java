
/*
 * @Author Chris Thatcher Jan 16. 2014
 */
public class MatrixTools {
	/*
	 * Static factory method to create the dimension X dimension Identity matrix. 
	 * 
	 * @param int the number of rows (alternatively, the number of columns) of the Identity matrix you want.
	 * @return A dimension X dimension identity matrix.
	 */
	public static Matrix makeIdentityMatrix(int dimension)
	{
		float[][] values = new float[dimension][dimension];
		
		for (int i = 0; i < dimension; i++)
		{
			for (int j = 0; j < dimension; j++)
			{
				values[i][j] = i==j ? 1 : 0;
			}
		}
		
		return new Matrix(values);
	}
	
	/*
	 * Static factory method to create a diagonally dominant and symmetric matrix of size "dimension"
	 * 
	 * @param float The lowest value that an element in the matrix may have.
	 * @param float The highest value that an element in the matrix may have.
	 * @param int the number of rows and columns in the matrix.
	 * @return A new randomized n by m matrix with values between lowerValueRange and upperValueRange that is 
	 * diagonally dominant and symmetric.
	 * 
	 */
	public static Matrix makeDiagonallyDominantMatrix(float lowerValueRange, float upperValueRange, int dimension)
	{
		float[][] values = new float[dimension][dimension];
		
		//make random matrix
		for (int i = 0; i < dimension; i++)
		{
			float rowSum = 0;
			
			for (int j = 0; j < dimension; j++)
			{
				values[i][j] = (float) ((Math.random() * (upperValueRange - lowerValueRange)) + lowerValueRange);
			}
		}
		
		//make them symmetric
		
		for (int i = 0; i < dimension; i++)
		{
			for (int j = 0; j < i; j++)
			{
				values[i][j] = values[j][i];
			}
		}
		
		//make it diagonally dominant
		
		for (int i = 0; i < dimension; i++)
		{
			float rowSum = 0;
			for (int j = 0; j < dimension; j++)
			{
				rowSum += Math.abs(values[i][j]);
			}
			values[i][i] = (float) (rowSum * 1.1 * (Math.random() > 0.5 ? 1 : -1));
		}
		
		return new Matrix(values);
	}
	
	/*
	 * Creates a new matrix that is the equivalent of augmenting Matrix B onto Matrix A. The original matrices remain unchanged.
	 * 
	 * @param Matrix The matrix that will be on the left side of the augmented matrix
	 * @param Matrix The matrix that will be on the right side of the augmented matrix
	 * 
	 * @return A new matrix that is equivalent to augmenting matrix A with matrix A.
	 */
	public static Matrix augmentMatrix(Matrix A, Matrix B)
	{
		float[][] aVal = A.getValues();
		float[][] bVal = B.getValues();
		
		int aRows = aVal.length;
		int bRows = bVal.length;
		
		if (aRows != bRows)
		{
			System.out.println("Cannot augment a matrix with a matrix with fewer rows\n" +
					"A rows = " + aRows + ", B rows = " + bRows);
		}
		
		int aCol = aVal[0].length;
		int bCol = bVal[0].length;
		
		int numRows = aRows;
		int numCol = aCol + bCol;
		
		
		float[][] values = new float[numRows][numCol];
		
		for (int i = 0; i < aRows; i++)
		{
			for (int j = 0; j < aCol; j++)
			{
				values[i][j] = aVal[i][j];
			}
		}
		
		for (int i = 0; i < bRows; i++)
		{
			for (int j = 0; j < bCol; j++)
			{
				values[i][j + aCol] = bVal[i][j];
			}
		}
		
		return new Matrix(values);
	}
	

	/*
	 * Static factory method to create a Matrix. Note that the matrix is build as float[rows][columns]
	 * 
	 * @param float[][] The values of the matrix. 
	 * @return A new matrix with data corresponding to values.
	 */
	public static Matrix makeMatrix(float[][] values)
	{
		return new Matrix(values);
	}
	
	/*
	 * Static factory method to create a randomized matrix. Elements in the matrix will be uniformly 
	 * distributed between lowerValueRange and upperValueRange inclusive.
	 * 
	 * @param float The lowest value that an element in the matrix may have.
	 * @param float The highest value that an element in the matrix may have.
	 * @param int The number of rows in the matrix.
	 * @param int the number of columns in the matrix.
	 * @return A new randomized n by m matrix with values between lowerValueRange and upperValueRange 
	 *
	 */
	public static Matrix makeRandomMatrix(float lowerValueRange, float upperValueRange, int rows, int columns)
	{
		float[][] tempValues = new float[rows][columns];
		
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < columns; j++)
			{
				tempValues[i][j] = (float) ((Math.random() * (upperValueRange - lowerValueRange)) + lowerValueRange);
			}
		}
		
		return new Matrix(tempValues);
	}
	
	/*
	 * Computes the product of two Matrices, A and B. Note that matrix 
	 * multiplication is not commutative, and thus the order of arguments 
	 * matters. In this case, the matrix on the left of the multiplication 
	 * should be passed in first, with the matrix on the right of the multiplication passed in second.
	 * 
	 * Note: This method works for both Matrices and Vectors, as Vectors are a 
	 * subclass of matrices. The method will return a Matrix, but the user can then call
	 * myReturnedMatrix.getVector() to get the corresponding vector. The getVector() call
	 * will throw an UnexpectedException if the return matrix is not a vector.
	 * 
	 * Note: If you multiply two vectors such that your result is a scalar, you will still get 
	 * a Matrix object back. You just need to call returnedMatrix.getScalar() if this is the case. 
	 * Note that that method just returns the float at [0][0] without checking, so the onus is 
	 * on you to verify you're doing it correctly.
	 * 
	 * @param Matrix The Matrix (or Vector) on the left of the multiplication
	 * @Param Matrix The Matrix (or Vector) on the right of the multiplication
	 * @return the product of AB.
	 */
	public static Matrix MatrixMultiplication(Matrix A, Matrix B)
	{	
		if (A.getColumns() != B.getRows())
		{
			System.out.println("Invalid Matrix Multiplcation:\nA.getColumns() = " + A.getColumns() + ", B.getRows() = " + B.getRows());
		}
		
		int resultingRows = A.getRows();
		int resultingColumns = B.getColumns();
		int innerCount = A.getColumns();
		
		
		float[][] result = new float[resultingRows][resultingColumns];
		
		float[][] aData = A.getValues();
		float[][] bData = B.getValues();
		
		for (int row = 0; row < resultingRows; row++)
		{
			for (int column = 0; column < resultingColumns; column++)
			{
				float value = 0;
				
				for (int count = 0; count < innerCount; count++)
				{
					value += aData[row][count] * bData[count][column];
				}
				
				result[row][column] = value;
			}
		}
		
		return new Matrix(result);
	}
	
	
	/*
	 * Uses Jacobi iteration to solve Ax = b. Given A and b, returns the vector x.
	 * 
	 * @param Matrix A, the matrix A.
	 * @param Vector b, the result of multiplying A and x
	 * 
	 * @return the Vector x.
	 */
	public static Vector solveForX(Matrix A, Vector b, float tolerance, int maxIter)
	{
		int count = 0;
		
		int rows = A.getRows();
		
		Vector v0 = MatrixTools.makeNonZeroVector(rows); 
		Vector v1 = v0;
		
		while (count < maxIter && (v1.equals(v0) || computeVectorError(v0, v1) > tolerance))
		{
			v0 = v1;
			
			Matrix diagonal = A.getDiagonal();
			Matrix remainder = Matrix.subtract(A, diagonal);
			
			Matrix inverseDiagonal = A.getInverseDiagonal();
			
			Matrix rx = MatrixTools.MatrixMultiplication(remainder, v0);
			Matrix bMinusRx = Matrix.subtract(b, rx);
			
			v1 = MatrixTools.MatrixMultiplication(inverseDiagonal, bMinusRx).getVector();
			count++;
		}
		
		//System.out.println("count = " + count);
		return v1;
	}
	
	
	/*
	 * Uses Jacobi iteration to solve Ax = b. Given A and b, returns the vector x.
	 * Uses GPU acceleration to, you know, accelerate stuff.
	 * 
	 * @param Matrix A, the matrix A.
	 * @param Vector b, the result of multiplying A and x
	 * 
	 * @return the Vector x.
	 */
	public static Vector solveForXGPUStyle(Matrix A, Vector b, float tolerance, int maxIter)
	{
		int count = 0;
		
		int rows = A.getRows();
		
		Vector v0 = MatrixTools.makeNonZeroVector(rows); 
		Vector v1 = v0;
		
		while (count < maxIter && (v1.equals(v0) || computeVectorError(v0, v1) > tolerance))
		{
			v0 = v1;
			
			Matrix diagonal = A.getDiagonal();
			Matrix remainder = Matrix.subtract(A, diagonal);
			
			Matrix inverseDiagonal = A.getInverseDiagonal();
			
			MatrixMultiplier k1 = new MatrixMultiplier(remainder, v0);
			k1.execute();
			Matrix rx = new Matrix(k1.getData());
			k1.dispose();
			
			Matrix bMinusRx = Matrix.subtract(b, rx);
			
			MatrixMultiplier k2 = new MatrixMultiplier(inverseDiagonal, bMinusRx);
			k2.execute();
			v1 = new Matrix(k2.getData()).getVector();
			k2.dispose();
			
			count++;
		}
		
		//System.out.println("count = " + count);
		return v1;
	}
	
	/*
	 * returns the percent error between two vectors.
	 */
	public static float computeVectorError(Vector oldVector, Vector newVector)
	{
		float differenceMag = (Matrix.subtract(newVector, oldVector)).getVector().getMagnitude();
		return (float) (Math.abs(differenceMag) / Math.abs(oldVector.getMagnitude()));
	}
	
	/*
	 * returns a nonzero vector of size "rows";
	 */
	public static Vector makeNonZeroVector(int rows)
	{
		float[] values = new float[rows];
		
		for (int i = 0; i < rows; i++)
		{
			values[i] = (float) (Math.random() + 0.01);
		}
		
		return new Vector(values, Orientation.COLUMN);
	}
}
