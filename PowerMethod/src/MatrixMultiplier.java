import com.amd.aparapi.Range;


public class MatrixMultiplier {
	private final Matrix m1;
	private final Matrix m2;
	private final Range range;
	private float[] linearResult;
	
	private MatrixMultiplierKernel k;
	
	public MatrixMultiplier(Matrix m1, Matrix m2)
	{
		this.m1 = m1;
		this.m2 = m2;
		int innerCount = m1.getColumns();
		int outerCount = m2.getColumns();
		
		this.linearResult = new float[m1.getRows() * m2.getColumns()];
		range = Range.create2D(m1.getRows(), m2.getColumns());
		
		k = new MatrixMultiplierKernel(getLinearData(m1.getValues()), getLinearData(m2.getValues()), linearResult, innerCount, outerCount);
	}
	
	private float[] getLinearData(float[][] data)
	{
		float[] linearData = new float[data.length * data[0].length];
		
		int count = 0;
		
		for (int i = 0; i < data.length; i++)
		{
			for (int j = 0; j < data[i].length; j++)
			{
				linearData[count++] = data[i][j];
			}
		}
		
		return linearData;
	}
	
	public float[][] getData() 
	{
		float[][] goodData = new float[m1.getRows()][m2.getColumns()];
		
		int counter = 0;
		for (int i = 0; i < m1.getRows(); i++)
		{
			for (int j = 0; j < m2.getColumns(); j++)
			{
				goodData[i][j] = linearResult[counter++];
			}
		}
		
		return goodData;
	}
	
	public void execute () {
		k.execute(range);
	}
	
	public void dispose() {
		k.dispose();
	}
}
