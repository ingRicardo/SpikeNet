/**
 * 
 */
package nnets;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;

import javax.imageio.ImageIO;

/**
 * @author TIJUANA
 *
 */
public class KMeansClusterV2 {

	/**
	 * @param args
	 */
	static int x1, y1, z1;
	static ArrayList<Integer> tmp = new ArrayList<Integer>();
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
	
		BufferedImage image = null;
		
		File input_file = new File("C:/Users/RICA/Documents/workspace-spring-tool-suite-4-4.0.0.RELEASE/JavaProjects/src/nnets/woman.jpg");
		
	
				
		try {
			image = ImageIO.read(input_file);
			
			x1 = image.getWidth();
			y1 = image.getHeight();
			z1 = 3;	
			
			System.out.println("image.getHeight() ->> "+ image.getHeight());
			System.out.println("image.getWidth() ->> "+ image.getWidth());
			
		//	image = new BufferedImage(x1, y1, BufferedImage.TYPE_INT_ARGB);
			
			
			double[][][] imag = new double[x1][y1][z1];
		    for (int i = 0; i < x1; i++) {
		        for (int j = 0; j < y1; j++) {
		        	for (int z = 0; z < z1; z++) {
		        		imag[i][j][z] = image.getRGB(i, j);
		        	//	System.out.println("imag[i][j][z] - >" + imag[i][j][z]);
		        	}
		        }
		    }
			
			
			
		//	double [][][] imag = new double[x1][y1][z1];
			
		    
		    
			double [][] imagx = new double[x1*y1][z1];
			int conta = 0;
			
			
			for (int i=0; i < x1; i++ ) {
				
				for (int j =0; j< y1; j++) {
					
					for (int k=0; k < z1; k++) {
						
						imagx[conta][k] = imag[i][j][k]; 
					//	System.out.println(imagx[conta][k]);
					}
					conta++;
				}
				
				
			}
			
			System.out.println("Done!");
			
			//Kmeans_var(imagx, 4, 500, 0.01);
			Kmeans_var(imagx, 4, 2, 0.01);
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// WRITE IMAGE 
      /*  try
        { 
            // Output file path 
            File output_file = new File("C:/projects/JavaProjects/src/nnets/chubby.jpg"); 
  
            // Writing to file taking type and path as 
            ImageIO.write(image, "jpg", output_file); 
  
            System.out.println("Writing complete."); 
        } 
        catch(IOException e) 
        { 
            System.out.println("Error: "+e); 
        } 
		
		*/
	}
	
	public static void Kmeans_var(double [][] xin, int nc, int maxepoca , double tol ) {
		
		double [][] imagy = new double[x1*y1][z1];
		
		//double dJ = tol+1;
		double dJ = 0;
		int epoca = 0;
		int lin = xin.length;
		int col = xin[0].length;
		
		double [][] input  =  xin;
		
		
		Random r = new Random();
		
		double [] ndx = new double[lin];
		
		for (int  i = 0; i< nc; i++) {
			ndx[i] = Math.round( (0 + (1 - nc) * r.nextDouble()) );
			System.out.println("ndx[ "+i+" ] -> "+ ndx[i]);
		}
			
		double [][] c = new double [nc][col];
		
		for (int n=0; n < nc ; n++) {
			
			for (int j = 0 ; j < col ; j++) {
				input[n][j] = ndx[n]; 
				c[n][j] = input[n][j];
				
				System.out.println("c[ "+n+" ][ "+j+" ] - > "+c[n][j]);
			}
			
			
		}
		 
		double [][] u = null;//new double[lin][nc];
		double [] dist = new double [nc];
		double sum =0;
		int um;
		double [] soma2 =new double [col];
		//double [] u = ;
		double [] dist1 = new double[nc];
		System.out.println("dJ -> "+ dJ + " tol "+ tol + " epoca "+ epoca + " maxepoca " + maxepoca);
		double soma1=0;
		double [] J= new double [maxepoca];
		while(dJ < tol && epoca< maxepoca) {
			System.out.println("inside while loop");
			u = new double[lin][nc];
			for(int vector =0; vector < lin; vector++) {
				
				for (int centro = 0; centro < nc; centro++) {
					//dist[centro] = (input[vector])
					sum = 0;
					for (int j =0; j < col ; j++) {
						
						sum = sum + ( input[vector][j] - Math.pow(c[centro][j], 2));
						System.out.println("sum -> "+ sum);
					}
					dist[centro] = sum;
				}
				
				 um=Find_min(dist);
					System.out.println(" um -> "+ um);
				u[vector][um] =1;
			}
			
			System.out.println(" u -> row- > "+ u.length +" lin "+ lin+ " col "+ u[0].length + " nc "+ nc);
			System.out.println("nc -> "+ nc);
			for (int centro = 0; centro < nc ; centro++) {
				double [] tmp = new double[lin];
				int [] ones = new int[lin];
				for (int p=0; p< lin; p++) {
					System.out.println("u[ "+lin+" ][ "+centro+" ]");
					tmp [p] = u[p][centro];

				}
				ones=find(tmp,1);
				System.out.println("<-ones->");
				print(ones);
				soma1 =0;
				
				soma2 =new double [col];
				double sumt=0;
				double sumV=0;
				for (int membro = 0; membro < ones.length; membro++) {
					//soma1 =soma1 +
					
				//System.out.println("ones[ "+membro+" ] - >"+ones[membro]);	
					soma1 = 0;
					sumt = 0;
					for (int x = 0; x < col; x++) {
						
						sumV =0;
					System.out.println(" input[ones[ "+membro+" ]][ "+x+" ] - > "+ input[ones[membro]][x]);	
					for (int y = 0; y< c[0].length; y++) {
						System.out.println("c[ "+centro+" ][ "+y+" ] -> "+c[centro][y]);
						
						sumt = Math.pow(input[ones[membro]][x]-c[centro][y], 2);
						soma1 = soma1 + sumt;
						System.out.println(" x-> "+ x);
					//	System.out.println("soma2[1][ "+x+" ] -> "+ soma2[1][x]);
						//soma2[1][x]= soma2[1][x] +input[ones[membro]][x];
						sumV = sumV + input[ones[membro]][x];
						
						
					}
					soma2[x] = sumV;
					System.out.println();
					}
				}
				
				///
				dist1[centro] = soma1;
				for (int y = 0; y< col; y++) {
					c[centro][y] = soma2[y]/ones.length;
				}
			}
			
			System.out.println("break");
			//break;
			epoca++;
			double sumds = 0;
			
			for (int x = 0; x < dist1.length; x++ ) {
				sumds = sumds + dist[x];
				
			}
			J[epoca] = sumds;
			
			if (epoca ==1) {
				dJ=J[epoca];
			}else {
				
				dJ = Math.abs(J[epoca-1]-J[epoca]);
			}
			
		}
		///////////////////////////
		double [] saida = new double[input.length];
		int [] clas = new int [lin];
		double [] clatmp =  new double[lin];
		for (int k=0; k< lin; k++) {
			clatmp =  new double[lin];
			for (int ck=0; ck< u[0].length; ck++) {
				clatmp[ck]= u[k][ck]; 
			}
			clas =findNotEqual(clatmp, 0);
		}
	}
	
	public static void print( int arr[]) {
		for (int d : arr) {
			System.out.print("   "+ d);
		}
	}
	public static int [] find(double [] vec, int val) {
	
		
		for (int j =0; j < vec.length; j++) {
			if (vec[j] == val)
				tmp.add(j);
			
		}
		int [] ones = new int[tmp.size()];
		for (int x = 0; x< tmp.size(); x++) {
			ones[x] = tmp.get(x);
		}
		
		return ones;	
	}
	
	public static int [] findNotEqual(double [] vec, int val) {
		
		for (int j =0; j < vec.length; j++) {
			if (vec[j] != val)
				tmp.add(j);
			
		}
		int [] ones = new int[tmp.size()];
		for (int x = 0; x< tmp.size(); x++) {
			ones[x] = tmp.get(x);
		}
		
		return ones;
		
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
	
	
public static int Find_min(double [] vec) {
		
		double ele = vec[0];
			double min =0;
			int indx =0;
			for (int x =1 ; x < vec.length; x++) {
				double vectmp = vec[x];
				
				if (ele < vectmp) {
					min = ele;
					indx = x;
				}else {
					min = vectmp;
					indx = x;
				}
									
				ele = min;
				
			}
				
	return indx;

	}

}
