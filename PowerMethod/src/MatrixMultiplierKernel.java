
public class MatrixMultiplierKernel extends com.amd.aparapi.Kernel {

	private float[] finalLinearData;
	private final float[] m1Data;
	private final float[] m2Data;
	private final int innerColumnCount;
	private final int outerColumnCount;
	
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

}
