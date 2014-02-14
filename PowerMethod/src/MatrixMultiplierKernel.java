
public class MatrixMultiplierKernel extends com.amd.aparapi.Kernel {

	private float[] finalLinearData;
	private float[] m1Data;
	private float[] m2Data;
	private int innerColumnCount;
	private int outerColumnCount;
	
	public MatrixMultiplierKernel(float[] m1LinearData, float[] m2LinearData, float[] finalLinearData, int innerColumnCount, int outerColumnCount)
	{
		this.m1Data = m1LinearData;
		this.m2Data = m2LinearData;
		this.finalLinearData = finalLinearData;
		this.innerColumnCount = innerColumnCount;
		this.outerColumnCount = outerColumnCount;
	}
	
	public float[] getData() {
		return this.finalLinearData;
	}
	
	@Override
	public void run() {
		int row = getGlobalId(0);
		int col = getGlobalId(1);
				
		for (int i = 0; i < innerColumnCount; i++)
		{
			finalLinearData[(row * outerColumnCount) + col] += m1Data[(row * innerColumnCount) + i] * m2Data[(i * outerColumnCount) + col];
		}
	}
	
	public void setM1LinearData(float[] m1)
	{
		this.m1Data = m1;
	}
	
	public void setM2LinearData(float[] m2)
	{
		this.m2Data = m2;
	}
	
	public void setFinalLinearData(float[] m3)
	{
		this.finalLinearData = m3;
	}

	public void setInnerColumnCount(int count)
	{
		this.innerColumnCount = count;
	}
	
	public void setOuterColumnCount(int count)
	{
		this.outerColumnCount = count;
	}
}
