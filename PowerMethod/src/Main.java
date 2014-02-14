import java.rmi.UnexpectedException;
import java.sql.SQLException;
import java.util.Arrays;

/*
 * @Author Chris Thatcher Jan 16. 2014
 */
public class Main {
	
	/*
	 * this main function tests out old school solutions vs. GPU accelerated solutions
	 */
	public static void main(String[] args) throws UnexpectedException, SQLException
	{
		Matrix A = MatrixTools.makeRandomMatrix(0, 1, 800, 800);
		Matrix B = MatrixTools.makeRandomMatrix(0, 1, 800, 800);
			
		
		System.out.println("Testing matrix multiplication: Two matrices of size 800 X 800");
		long start = System.currentTimeMillis();
		Matrix C = MatrixTools.MatrixMultiplication(A, B);
		long end = System.currentTimeMillis() - start;
		System.out.println("Old School took " + end);
		//C.print();
		
		start = System.currentTimeMillis();
		MatrixMultiplier k = new MatrixMultiplier(A, B);
		k.execute();
		
		Matrix D = new Matrix(k.getData());
		end = System.currentTimeMillis() - start;
		
		System.out.println("GPU Accelerated took: " + end);
		//D.print();
		
		Matrix E = Matrix.subtract(C, D);
		
		//E.print();
		double d   = 0;
		for (float[] f : E.getValues())
		{
			for (float fl : f)
			{
				d += fl;
			}
		}
		
		System.out.println("Difference is " + d);
		
		System.out.println();
		
		System.out.println("Testing: Find largest eigenvalue on 5000 X 5000 matrix");
		A = MatrixTools.makeDiagonallyDominantMatrix(0, 10, 5000);
		
		start = System.currentTimeMillis();
			float old = A.getGreatestEigenValueSmart(0.0005f, 50000);
		end = System.currentTimeMillis() - start;
		
		System.out.println("Old greatest Eigenvalue took: " + end + " and was " + old);
		
		start = System.currentTimeMillis();
			float nuevo = A.getGreatestEigenValueSmartGPUStyle(0.0005f, 50000);
		end = System.currentTimeMillis() - start;
		
		System.out.println("GPU Greatest Eigenvalue took: " + end + " and was " + nuevo);
		
		System.out.println("Error is: " + MathUtils.computeRelativeError(old, nuevo));
		
		System.out.println();
		
		System.out.println("Testing: Find smallest eigenvalue on 2500 X 2500 matrix");
		A = MatrixTools.makeDiagonallyDominantMatrix(0, 10, 2500);
		
		start = System.currentTimeMillis();
			old = A.getSmallestEigenValueUsingJacobi(0.0005f, 50000);
		end = System.currentTimeMillis() - start;
		
		System.out.println("Old smallest Eigenvalue took: " + end + " and was " + old);
		
		start = System.currentTimeMillis();
			nuevo = A.getSmallestEigenValueUsingJacobiGPUStyle(0.0005f, 50000);
		end = System.currentTimeMillis() - start;
		
		System.out.println("GPU smallest Eigenvalue took: " + end + " and was " + nuevo);
		
		System.out.println("Error is: " + MathUtils.computeRelativeError(old, nuevo));

		
		System.out.println();
		System.out.println("Testing: Find condition number on a 2500 X 2500 matrix");
		A = MatrixTools.makeDiagonallyDominantMatrix(0, 10, 2500);
		
		start = System.currentTimeMillis();
			old = A.getConditionNumber();
		end = System.currentTimeMillis() - start;
		
		System.out.println("Old condition number took: " + end + " and was " + old);
		
		start = System.currentTimeMillis();
			nuevo = A.getConditionNumberGPUStyle();
		end = System.currentTimeMillis() - start;
		
		System.out.println("GPU condition number took: " + end + " and was " + nuevo);
		
		System.out.println("Error is: " + MathUtils.computeRelativeError(old, nuevo));
		
	}
	
	/*
	 * Computes the condition number of "count" different randomly generated symmetric diagonally dominant matrices.
	 * Displays the median condition number and average condition number. 
	 * Matrices are of size 50, with normal values ranging between -10 and 10, with diagonal values guaranteed to be 
	 * greater in magnitude than the sum of the magnitudes of each value in the row.
	 * 
	 * 
	 */
	private static void runConditionNumberTests(int count) 
	{
		float[] results = new float[count];
		for (int i = 0; i < count; i++)
		{
			Matrix A = MatrixTools.makeDiagonallyDominantMatrix(-10, 10, 50);
			results[i] = A.getConditionNumber();
		}
		
		Arrays.sort(results);
		System.out.println("Median Condition Number = " + results[results.length / 2]);
		float sum = 0;
		for (int i = 0; i < results.length; i++)
		{
			sum += results[i];
		}
		System.out.println("Average Condition Number = " + sum / results.length);
	}
	
	/*
	 * Creates "count" Diagonally dominant symmetric matrices of size 50. 
	 * Computes the smallest Eigenvalue of that matrix using Jacobi Iteration as well as using the real inverse of the matrix. 
	 * Displays the Median difference between those two values.
	 */
	private static void runJacobiIterationTests(int count) 
	{
		//Run count tests, and check the median difference between the smallestEigenValue using inverse and smallestEigenValue using Jacobi
		
		float[] results = new float[count];
		
		for (int i = 0; i < count; i++)
		{
			Matrix A = MatrixTools.makeDiagonallyDominantMatrix(-10, 10, 50);
			results[i] = A.getSmallestEigenValue(0.00005f, 50000) - A.getSmallestEigenValueUsingJacobi(0.00005f, 50000);
		}
		Arrays.sort(results);
		System.out.println("Median difference between jacobi smallest eigenvalue and real inverse smallest eigenvalue:");
		System.out.println(results[results.length / 2]);
	}
	

	/*
	 * This method compares the runtime of the smart eigenvalue finding algorithm vs. the naive eigenvalue finding algorithm vs. the Inverse Power Method.
	 */
	private static void runPowerMethodTests() throws UnexpectedException, SQLException
	{
		for (int i = 20; i < 401; i+=5)
		{
			for (int j = 0; j < 20; j++)
			{
				Matrix temp = MatrixTools.makeRandomMatrix(-100, 100, i, i);
				
				//do naive test
				long start = System.currentTimeMillis();
				float result = temp.getGreatestEigenvalueNaive(0.0005f, 10000);
				long totalDur = System.currentTimeMillis() - start;
				
//				System.out.println("Naive Result = " + result + ", in " + totalDur + " milliseconds");
				
				//do smart test
				start = System.currentTimeMillis();
				result = temp.getGreatestEigenValueSmart(0.0005f, 10000);
				totalDur = System.currentTimeMillis() - start;
				
//				System.out.println("Smart Result = " + result + ", in " + totalDur + " milliseconds");
				
				//Do smallest test
				start = System.currentTimeMillis();
				result = temp.getSmallestEigenValue(0.0005f, 10000);
				totalDur = System.currentTimeMillis() - start;
				
				System.out.println("Inverse Result = " + result + ", in " + totalDur + " milliseconds");
				
//				System.out.println("\n\n");

			}
		}
	}


 
}
