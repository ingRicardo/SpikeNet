/**
 * 
 */
package nnets;

/**
 * @author RICA
 *
 */
public class SpikeTest {

	/**
	 * 
	 */
	public SpikeTest() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		double ein[][] = {{21.0,9.0,13.0,4.0,6.0,8.0,0.0,5.0,20.0},
					   {8.0,7.0,14.0,4.0,15.0,5.0,3.0,1.0,2.0}};
		
		double maxd = 0.2;
		double rho = 0.5;
		boolean show = false;
		double range =0.0;
		//double einnorma[][] = new double[ein[0].length][ein[0].length];
		
		System.out.println("print Data");
		for (double[] is : ein) {
			System.out.println(" ");
			print(is);
		}
		
		
		int r=0;
		for (double[] row : ein) {
			
			range = MaxInArray(row) - MinInArray(row);
			double[] rown = new double[ein[0].length];
			for (int x =0; x < row.length; x++) {
				rown[x] = (row[x] - MinInArray(row))/range;	
			}
			ein[r] = rown;
			r++;
		}
	
		System.out.println(" ");
		System.out.println("Normalized Data");
		for (double[] is : ein) {
			System.out.println(" ");
			print(is);
		}
		
		}
	
	public static void print(double[] vec ) {
		
		for (int i=0; i < vec.length; i++) {
			
			System.out.println("val vec ->> "+ vec[i]);
		}
	}
	
	public static double MaxInArray(double vec[]) {
		
		double may =0;
		double run = 0;
		double ini = vec[0];
		
			for (int i =1; i < vec.length; i++) {
				run = vec[i];
				
				if (ini > run )
					may = ini;
				else
					ini =run;
				
				if (ini == run)
					may =ini;	
			}
		
		return may;
	}
	
public static double MinInArray(double vec[]) {
		
		double min =0;
		double run = 0;
		double ini = vec[0];
		
			for (int i =1; i < vec.length; i++) {
				run = vec[i];
			
			if (ini != run)
				if (ini > run )
					ini =run;
				else
					min = ini;
			
			if (ini == run)
				min =ini;		
				
			}
		
		return min;
	}
	
	public static void kodieren_ln(int ein[][], double maxd, double rho, boolean plot) {
		
		
		
	}
	
}
