/**
 * 
 */
package nnets;

import java.util.Random;

/**
 * @author RICA
 *
 */
public class Clusterv2 {



	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		test();
	}
	
	public static void print( double arr[]) {
		for (double d : arr) {
			System.out.print("   "+ d);
		}
	}
	
	public static double max(double [] vec) {
			double ele = vec[0];
				double max =0;
				for (int x =1 ; x < vec.length; x++) {
					double vectmp = vec[x];
					
					if (ele > vectmp) {
						max = ele;
					}else 
						max = vectmp;
					
					ele = max;
				}
					
		return max;
	}
	
	public static double min (double [] vec) {
		
		double ele = vec[0];
			double min =0;
			for (int x =1 ; x < vec.length; x++) {
				double vectmp = vec[x];
				
				if (ele < vectmp) {
					min = ele;
				}else 
					min = vectmp;
				
				ele = min;
			}
				
	return min;

	}
	public static double [][] kodieren_rf(double [][] ein, int rf[], 
			int maxd, double [] sigma, double f_cut, double rho, int typ, int plt) {
		
		int t_neur = 0;
		for (int i =0; i < rf.length; i++) {
			t_neur = t_neur + rf[i]; 
			
		}
		
		int t_cj =  ein.length;
		int t_in =  ein[0].length;
		
		double [][]  aus = new double [t_cj][ t_neur]; 
		
		
		int n=0;
		
		double [] ct = new double[rf.length];
		if (rf.length <3)
			System.out.println("receptive fields must be greater than 3");
		
		System.out.println(" ");
		double range =0;
		int i = 0;
		for (double [] d : ein) {
			range =0;
			range =max(d) -min(d);  
			System.out.println(" max "+ max(d)+" min "+ min(d)+" range "+ range);
			for (int j = 0; j < d.length; j++) {				
				ein[i][j]=	(ein[i][j] - min(d))/range;
			}
			  
			i++;
		}
		
	//	double [] ct = new double[rf.length];
		
		if (typ ==1) {
			for (int j =0 ; j < rf.length; j++ ) {
				ct[j] = (j-0.5) * (1/rf[j]);
			}
		}else if (typ == 0) {
				for (int j =0 ; j < rf.length; j++ ) {
					ct[j] = (j-1.5) * (1/rf[j]-2);
				}
				
			}else
				System.out.println("typ must be 0 / 1");
			
		
		
		
		System.out.println("\nPRINT EIN NORMALIZED \n");
		for (double [] d : ein) {
			System.out.println(" ");
			print(d);
			
			
		}
		
		double aux = 0; 
		
		for (int j =0; j< t_cj; j++) {
			
			for (int k=0; k < rf.length; k++) {
				aux = maxd-maxd* Math.exp(-1*(Math.pow((ein[j][k]-ct[k]), 2)
						/ Math.pow(2*sigma[k], 2)));
				aus[j][k] = (aux*(1/rho))*rho;
				
				if (aus[j][k] > maxd*f_cut)
					aus[j][k] = -1;
				
				//System.out.println( ein[j][k]);
				//System.out.println( aus[j][k]);
			}
		}
		
	
		
		
		return aus;
	}
	public static void test() {
		int out_neu=5;
		double rho = 0.1;
		int conj = 5;
		int tau = 3;
		int [] rf = {1, 1, 1};
		double [] sigma = {1/(1.5*rf[0]-2),1/(1.5*rf[1]-2), 1/(1.5*rf[0]-2)};
		double f_cut =0.9;
		int maxd = 2;
		
		///
		int [] d = {0,1,2};
		int w_max = 1, w_min=0;
		int max_epoch = 2;
		int t_max =5;
		double eta = 0.25;
		double beta = 0.1;
		double nu = 5.0;
		double dx = 2.3;
		////
		
		double kapa = 1 - (Math.pow(nu, 2)) / (2* Math.log(beta/(beta+1)) ); 
		
		int in_neu = 0;
		for (int i =0 ; i < rf.length; i++) {
			in_neu = in_neu + rf[i];
			
		}
		
		System.out.println("in_neu -> "+ in_neu);
		int ssin = d.length;
		double [][][] w= new double[in_neu][out_neu][ssin];
		
	
		double [][] ein = {{2.3,6, 8.5},
						   {9.3,3.3, 4.2},
						   {15,3.5,7.2},
						   {22.1,6.6, 5.7},
						   {1.3,8.4,1.6}};
		double ssinv=0;
		Random r = new Random();
		for (int  i =0; i <in_neu; i++)
			for( int j=0; j < out_neu; j++) {
				 ssinv = 0;
				for (int s =0; s < ssin; s++) {
					 ssinv=	(0 + (1 - 0) * r.nextDouble());									
					 w[i][j][s] = (w_min+ssinv) * (w_max-w_min); 
				}
					
			}
		
		System.out.println("PRINT W \n");
		
		for (double[][] es : w) {
			System.out.println(" ");
			for (double[] es2 : es) {
				print(es2);
			}
			
		}
		
		// delays
				System.out.println("\n delays \n");
		
		double [][] aus = kodieren_rf(ein, rf, maxd, sigma, f_cut, rho, 0, 0);
		
		for (double [] d1 : aus) {
			System.out.println(" ");
			print(d1);
		}
			
	}

}
