import java.rmi.UnexpectedException;

/*
 * @Author Chris Thatcher Jan 16. 2014
 */
public class Matrix {
	private int rows;
	private int columns;
	protected float[][] values;
	private boolean isScalar = false;

	/*
	 * Class Constructor
	 * @param float[][] the elements of the matrix.
	 */
	public Matrix(float[][] values)
	{
		this.values = values;
		this.rows = values.length;
		this.columns = values[0].length;
		
		if (rows == 1 && columns == 1)
		{
			isScalar = true;
		}
	}
	
	/*
	 * Constructor to create a backing matrix from a vector. Given a vector an orientation, 
	 * we create a Matrix that either has one row of n elements or n rows of 1 elements. 
	 * 
	 * @param float[] the elements in the vector
	 * @param Orientation the orientation of the vector(ROW or COLUMN)
	 * 
	 * @return a new Matrix that represents the given Vector.
	 */
	public Matrix(float[] column, Orientation orientation)
	{
		float[][] values;
		
		if (orientation == Orientation.ROW)
		{
			values = new float[1][column.length];
			for (int i = 0; i < column.length; i++)
			{
				values[0][i] = column[i];
			}
		}
		else
		{
			values = new float[column.length][1];
			for (int i = 0; i < column.length; i++)
			{
				values[i][0] = column[i];
			}
		}
		
		this.values = values;
		this.rows = values.length;
		this.columns = values[0].length;
		
		if (rows == 1 && columns == 1)
		{
			isScalar = true;
		}
	}
	
	/*
	 * If this Matrix represents a vector, returns a new vector that is equivalent to the matrix. 
	 * If this matrix does not represent a vector (that is, it it not one dimensional) it will throw an 
	 * UnexpectedException.
	 * 
	 * @return New Vector that is equivalent to the Vector represented by this Matrix.
	 */
	public Vector getVector() {
		if (rows == 1)
		{
			return new Vector(values[0], Orientation.ROW);
		}
		else if (columns == 1)
		{
			float[] column = new float[values.length];
			for (int i = 0; i < values.length; i++)
			{
				column[i] = values[i][0];
			}
			
			return new Vector(column, Orientation.COLUMN);
		}
		else
		{
			System.out.println("This Matrix is NOT a vector");
			return null;
		}
	}
	
	/*
	 * Returns the element of the matrix at 0,0. This should only be used when this Matrix is the result 
	 * of a matrix multiplication of two vectors resulting in a scalar. Note that there is no checking to make sure this is okay.
	 * 
	 * @return float the element at 0,0 of the matrix. If this matrix is a scalar, then the value of that scalar.
	 */
	public float getScalar() {
		return values[0][0];
	}
	
	/*
	 * Returns an approximation of the condition number of the calling matrix.
	 * It is found by approximating the largest and smallest eigenvalues. 
	 * Note: this only works for symmetric diagonally dominant matrices.
	 */
	public float getConditionNumber() {
		return Math.abs( getGreatestEigenValueSmart(0.0005f, 5000) / getSmallestEigenValueUsingJacobi(0.0005f,  5000));
	}
	
	/*
	 * Returns an approximation of the condition number of the calling matrix.
	 * It is found by approximating the largest and smallest eigenvalues. 
	 * Note: this only works for symmetric diagonally dominant matrices.
	 */
	public float getConditionNumberGPUStyle() {
		return Math.abs( getGreatestEigenValueSmartGPUStyle(0.0005f, 5000) / getSmallestEigenValueUsingJacobiGPUStyle(0.0005f,  5000));
	}
	
	/*
	 * Approximates the greatest Eigenvalue of the calling Matrix. It is approximated using the power method. 
	 * This method is naive in the fact that it computes two matrix multiplications each iteration instead of just one. 
	 * 
	 * @param float the tolerance of our approximation
	 * @param int the maximum number of iterations to run
	 * 
	 * @return An approximation of the greatest Eigenvalue of the calling matrix as a float.
	 */
	public float getGreatestEigenvalueNaive(float tolerance, int maxIterations) throws UnexpectedException 
	{
		int count = 0;
		float error = tolerance * 10;
		float lambda0 = Float.MAX_VALUE;
		
		Matrix A = this;
		float[] guessVectorData = new float[rows];
		for (int i = 0; i < rows; i++)
		{
			if (i % 2 == 0)
			{
				guessVectorData[i] = 1;
			}
			else
			{
				guessVectorData[i] = 2;
			}
		}
		Vector z0 = Vector.makeVector(guessVectorData, Orientation.COLUMN);
		
		
		while (error > tolerance && count < maxIterations)
		{		
			Vector z1 = MatrixTools.MatrixMultiplication(A, z0).getVector();
			
			Vector y = z1.getNormalizedVector();
			
			float lambda1 = MatrixTools.MatrixMultiplication(y.transpose(), MatrixTools.MatrixMultiplication(A, y)).getScalar();
			
			error = MathUtils.computeRelativeError(lambda0, lambda1);
			z0 = y;
			lambda0 = lambda1;
			count++;
		}
		
		return lambda0;
	}
	
	/*
	 * Approximates the greatest Eigenvalue of the calling Matrix. It is approximated using the power method. 
	 * This method is smarter in the fact that it only computes one matrix multiplication per iteration. 
	 * 
	 * @param float the tolerance of our approximation
	 * @param int the maximum number of iterations to run
	 * 
	 * @return An approximation of the greatest Eigenvalue of the calling matrix as a float.
	 */
	public float getGreatestEigenValueSmart(float tolerance, int maxIterations)
	{
		int count = 0;
		float error = tolerance * 10;
		float lambda0 = Float.MAX_VALUE;
		float lambda1;
		
		Matrix A = this;
		float[] guessVectorData = new float[columns];
		for (int i = 0; i < columns; i++)
		{
			if (i % 2 == 0)
			{
				guessVectorData[i] = 1;
			}
			else
			{
				guessVectorData[i] = 2;
			}
		}
		Vector z0 = Vector.makeVector(guessVectorData, Orientation.COLUMN);
		
		Vector z1 = MatrixTools.MatrixMultiplication(A, z0).getVector();
		
		while (error > tolerance && count < maxIterations)
		{
			Vector y = z1.getNormalizedVector();
			Vector w = MatrixTools.MatrixMultiplication(A, y).getVector();
			
			lambda1 = MatrixTools.MatrixMultiplication(y.transpose(), w).getScalar();
			
			error = MathUtils.computeRelativeError(lambda0, lambda1);
			
			z1 = w;
			lambda0 = lambda1;
			count++;
		}
		
		return lambda0;
	}
	
	/*
	 * Approximates the greatest Eigenvalue of the calling Matrix. Does Matrix multiplication using the GPU. 
	 * It is approximated using the power method. 
	 * This method is smarter in the fact that it only computes one matrix multiplication per iteration. 
	 * 
	 * @param float the tolerance of our approximation
	 * @param int the maximum number of iterations to run
	 * 
	 * @return An approximation of the greatest Eigenvalue of the calling matrix as a float.
	 */
	public float getGreatestEigenValueSmartGPUStyle(float tolerance, int maxIterations)
	{
		int count = 0;
		float error = tolerance * 10;
		float lambda0 = Float.MAX_VALUE;
		float lambda1;
		
		Matrix A = this;
		float[] guessVectorData = new float[rows];
		for (int i = 0; i < rows; i++)
		{
			if (i % 2 == 0)
			{
				guessVectorData[i] = 1;
			}
			else
			{
				guessVectorData[i] = 2;
			}
		}
		Vector z0 = Vector.makeVector(guessVectorData, Orientation.COLUMN);
		
		MatrixMultiplier k = new MatrixMultiplier(A, z0);
		k.execute();
		Vector z1 = new Matrix(k.getData()).getVector();
		
		while (error > tolerance && count < maxIterations)
		{
			Vector y = z1.getNormalizedVector();
			Vector w = MatrixTools.MatrixMultiplication(A, y).getVector();
			
			k = new MatrixMultiplier(y.transpose(), w);
			k.execute();
			lambda1 = new Matrix(k.getData()).getScalar();
			
			error = MathUtils.computeRelativeError(lambda0, lambda1);
			
			z1 = w;
			lambda0 = lambda1;
			count++;
		}
		
		return lambda0;
	}
	
	/*
	 * Approximates the smallest Eigenvalue of the calling matrix. Uses the inverse power method
	 * to compute the smallest eigenvalue. 
	 * 
	 * @param float the tolerance of our approximation
	 * @param int the maximum number of iterations to run
	 * 
	 * @return An approximation of the smallest Eigenvalue of the calling matrix as a float.
	 */
	public float getSmallestEigenValue(float tolerance, int maxIterations)
	{
		Matrix inverse = this.getInverse();
		
		float eigenValue = inverse.getGreatestEigenValueSmart(tolerance, maxIterations);
		
		return (float) Math.pow(eigenValue, -1);
	}
	
	
	/*
	 * Approximates the smallest Eigenvalue of the calling matrix. Uses a modified inverse power method
	 * to compute the smallest Eigenvalue. Instead of computing the inverse of the Matrix A, this method
	 * uses Jacobi iteration to approximate it.
	 * 
	 * @param float the tolerance of our approximation
	 * @param int the maximum number of iterations to run
	 * 
	 * @return An approximation of the smallest Eigenvalue of the calling matrix as a float.
	 */
	public float getSmallestEigenValueUsingJacobi(float tolerance, int maxIterations)
	{
		int count = 0;
		float error = tolerance * 10;
		float lambda0 = Float.MAX_VALUE;
		float lambda1;
		
		Matrix A = this;
		Vector z0 = MatrixTools.makeNonZeroVector(this.getRows());

		Vector z1 = MatrixTools.solveForX(A, z0, tolerance, maxIterations);
		
		while (error > tolerance && count < maxIterations)
		{
			Vector y = z1.getNormalizedVector();
			Vector w = MatrixTools.solveForX(A, y, tolerance, maxIterations);
			
			lambda1 = MatrixTools.MatrixMultiplication(y.transpose(), w).getScalar();
			
			error = MathUtils.computeRelativeError(lambda0, lambda1);
			
			z1 = w;
			lambda0 = lambda1;
			count++;
		}
		
		return 1/lambda0;
	}

	
	/*
	 * Approximates the smallest Eigenvalue of the calling matrix. Uses GPU Acceleration. Uses a modified inverse power method
	 * to compute the smallest Eigenvalue. Instead of computing the inverse of the Matrix A, this method
	 * uses Jacobi iteration to approximate it.
	 * 
	 * @param float the tolerance of our approximation
	 * @param int the maximum number of iterations to run
	 * 
	 * @return An approximation of the smallest Eigenvalue of the calling matrix as a float.
	 */
	public float getSmallestEigenValueUsingJacobiGPUStyle(float tolerance, int maxIterations)
	{
		int count = 0;
		float error = tolerance * 10;
		float lambda0 = Float.MAX_VALUE;
		float lambda1;
		
		Matrix A = this;
		Vector z0 = MatrixTools.makeNonZeroVector(this.getRows());

		Vector z1 = MatrixTools.solveForXGPUStyle(A, z0, tolerance, maxIterations);
		
		while (error > tolerance && count < maxIterations)
		{
			Vector y = z1.getNormalizedVector();
			Vector w = MatrixTools.solveForXGPUStyle(A, y, tolerance, maxIterations);
			
			MatrixMultiplier k = new MatrixMultiplier(y.transpose(), w);
			k.execute();
			lambda1 = new Matrix(k.getData()).getScalar();
			
			error = MathUtils.computeRelativeError(lambda0, lambda1);
			
			z1 = w;
			lambda0 = lambda1;
			count++;
		}
		
		return 1/lambda0;
	}
	
	/*
	 * Returns a new matrix that is equivalent to the inverse of the calling matrix, if it exists. If it does not exist an exception is thrown.
	 * 
	 * @return A new matrix that is equivalent to the inverse of the calling matrix. 
	 */
	public Matrix getInverse()
	{
		if (this.rows != this.columns)
		{
			System.out.println("Exception: Cannot compute the inverse of a non-square matrix");
		}
		
		Matrix temp = MatrixTools.augmentMatrix(this, MatrixTools.makeIdentityMatrix(this.getRows()));
		temp = temp.executeGaussJordanElimination();
		Matrix inverse = Matrix.getSubMatrix(temp, 0, this.getColumns());
		return inverse;
	}
	
	/*
	 * Returns a new matrix that is equivalent to performing Gauss Jordan elimination on the calling matrix. 
	 * 
	 * @return a new matrix that is equivalent to performing Gauss Jordan elimination on the calling matrix.
	 */
	private Matrix executeGaussJordanElimination()
	{
		Matrix temp = new Matrix(getValues());
		
		for (int i = 0; i < rows; i++)
		{
			if (temp.getValues()[i][i] == 0)
			{
				temp.switchRowsForNonZeroEntryAt(i, i);
			}
			
			for (int j = i + 1; j < rows; j++)
			{
				temp.makeThisEntryZero(j, i, i);
			}
		}
		
		//So now temp should be in triangle form. So now we need to back substitute and pwn this thing.
		
		for (int i = rows -1; i >= 0; i--)
		{
			for (int j = i - 1; j >= 0; j--)
			{
				temp.makeThisEntryZero(j, i, i);
			}
		}
		
		//Sweet action. So, now it's rref. We just need to make them all 1's.
		
		for (int i = 0; i < rows; i++)
		{
			temp.scaleThisRowToOne(i, i);
		}
		
		return temp;
	}
	
	/*
	 * Switches the current row (given by i) with the first row below it that has a non-zero entry in position j.
	 */
	private void switchRowsForNonZeroEntryAt(int locationi, int locationj)
	{
		int originalToSwap = locationi;
		int newToSwap = -1;
		
		for (int i = locationi + 1; i < this.rows; i++)
		{
			if (values[i][locationj] == 0)
			{
				newToSwap = i;
				break;
			}
		}
		
		if (newToSwap == -1)
		{
			System.out.println("Cannot invert a rank-deficient matrix!");
		}
		
		float[] tempRow = values[originalToSwap];
		values[originalToSwap] = values[newToSwap];
		values[newToSwap] = tempRow;
	}
	
	/*
	 * Performs row-reduction on the target row using the original row such that the element specified by values[origRow][columnToKill] becomes 0.
	 * 
	 * @param int TargetRow, or the row we need to reduce
	 * @param int origRow, or the row we will use to reduce it
	 * @param int columnToKill, or the column that we want to become 0
	 */
	public void makeThisEntryZero(int targetRow, int origRow, int columnToKill)
	{
		float origValue = values[origRow][columnToKill];
		float valueToKill = values[targetRow][columnToKill];
		
		float coefficient = valueToKill / origValue;
		
		for (int i = 0; i < values[targetRow].length; i++)
		{
			values[targetRow][i] = values[targetRow][i] - (values[origRow][i] * coefficient);
		}
	}
	
	/*
	 * Performs a row-reduction by dividing the targetRow by a scalar such that the targetEntry becomes 1.
	 * 
	 * @param int targetRow, or the row we want to reduce
	 * @param int targetEntry, or the column of the element we want to become 1
	 */
	public void scaleThisRowToOne(int targetRow, int targetEntry)
	{
		float origValue = values[targetRow][targetEntry];
		
		for (int i = 0; i < this.columns; i++)
		{
			values[targetRow][i] = values[targetRow][i] / origValue;
		}
	}

	/*
	 * Static method that returns the submatrix of a given matrix. This is useful for un-augmenting an augmented matrix.
	 * 
	 * @param Matrix orig, or the original matrix that we want to take a part from
	 * @param int startRow, or the first row we want included in the submatrix
	 * @param int startCol, or the first column we want included in the submatrix.
	 * 
	 * @return A new matrix that is equivalent to the submatrix defined by startRow to endRow and startCol to endCol.
	 */
	public static Matrix getSubMatrix(Matrix orig, int startRow, int startCol)
	{
		float[][] oldValues = orig.getValues();
		int totalRows = orig.getRows();
		int totalColumns = orig.getColumns();
		
		int newRowCount = totalRows - startRow;
		int newColCount = totalColumns - startCol;
		
		float[][] newValues = new float[totalRows - startRow][totalColumns - startCol];
		
		for (int i = 0; i < newRowCount; i++)
		{
			for (int j = 0; j < newColCount; j++)
			{
				newValues[i][j] = oldValues[i + startRow][j + startCol];
			}
		}	
		
		return new Matrix(newValues);
	}
	
	/*
	 * Returns a new matrix that is empty except for the diagonal. The diagonal is the same
	 * as the diagonal of the calling matrix.
	 * 
	 * @return An empty matrix with the diagonal of the calling matrix.
	 */
	public Matrix getDiagonal() {
		float[][] values = new float[this.rows][this.columns];
		
		for (int i = 0; i < this.rows; i++)
		{
				values[i][i] = this.getValues()[i][i];
		}
		
		return new Matrix(values);
	}

	/*
	 * Returns a new matrix with zeros everywhere except the diagonal, which is the inverse of the calling matrix's diagonal.
	 */
	public Matrix getInverseDiagonal() {
		float[][] values = new float[this.rows][this.columns];
		
		for (int i = 0; i < this.rows; i++)
		{
				values[i][i] = (float) Math.pow(this.getValues()[i][i], -1);
		}
		
		return new Matrix(values);
	}
	
	/*
	 * Returns a new matrix that is the result of subtracting matrix B from matrix A.
	 */
	public static Matrix subtract(Matrix A, Matrix B)
	{
		int rows = A.getRows();
		int cols = A.getColumns();
		
		float[][] AValues = A.getValues();
		float[][] BValues = B.getValues();
		
		float[][] newValues = new float[rows][cols];
		
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				newValues[i][j] = AValues[i][j] - BValues[i][j];
			}
		}
		
		return new Matrix(newValues);
	}
	
	/*
	 * Convenience method to display a matrix}
	 */
	public void print()
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < columns; j++)
			{
				System.out.print(values[i][j] + ", ");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	public boolean isScalar() {
		return isScalar;
	}
	public int getRows() {
		return rows;
	}
	public void setRows(int rows) {
		this.rows = rows;
	}
	public int getColumns() {
		return columns;
	}
	public void setColumns(int columns) {
		this.columns = columns;
	}
	public float[][] getValues() {
		return values;
	}
	public void setValues(float[][] values) {
		this.values = values;
	}
	
}
