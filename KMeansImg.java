/**
 * 
 */
package ai;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Random;

import javax.imageio.ImageIO;

/**
 * @author TIJUANA
 *
 */
public class KMeansImg {

	/**
	 * 
	 */
	private int[][] myPixels;
	
	 /**
     * Constructor takes in the file that contains the image
     * and converts it to a BufferedImage. This image is stored
     * in an instance variable in the form of an integer array.
     *
     *         DO NOT MODIFY THIS METHOD
     * 
     * @param file, file that contains the image
     */
	
	
	public KMeansImg(File file) {
		// TODO Auto-generated constructor stub
        BufferedImage img;
        try {
            img = ImageIO.read(file);
        } catch (Exception e) {
            System.out.println("Incorrect File ");
            return;
        }
        if (img.getType() == BufferedImage.TYPE_INT_RGB) {
            myPixels = imageToPixels(img);
        } else {
            BufferedImage tmpImage = new BufferedImage(img.getWidth(null), img.getHeight(null), BufferedImage.TYPE_INT_RGB);
            tmpImage.createGraphics().drawImage(img, 0, 0, null);
            myPixels = imageToPixels(tmpImage);
        }   
	}
	
    /**
     * Method to return the original image
     *
     *         DO NOT MODIFY THIS METHOD
     * 
     * @return image, the original image
     */
    public int[][] getPixels() {
        return myPixels;
    }
	
	/**
     * Converts a BufferedImage to a two-dimensional array of image data.
     *
     *         DO NOT MODIFY THIS METHOD
     *
     * @param image the source BufferedImage
     *
     * @return  a two-dimensional array of image data representing the
     *          source image
     *
     * @throws  IllegalArgumentException if the source image is
     *          ill-defined.
     */
    public static int[][] imageToPixels(BufferedImage image) throws IllegalArgumentException {
        if (image == null) {
            throw new IllegalArgumentException();
        }
        
        int width = image.getWidth();
        int height = image.getHeight();
        int[][] pixels = new int[height][width];
        for (int row = 0; row < height; row++) {
            image.getRGB(0, row, width, 1, pixels[row], 0, width);
        }
        return pixels;
    }
	/**http://www.cs.utexas.edu/users/mckinley/305j/lectures/BlueJExamples/24-color/Transformations.java
	 * @param args
	 * @throws IOException 
	 * https://stackoverflow.com/questions/21576096/convert-a-2d-array-of-doubles-to-a-bufferedimage
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		File file = new File("/home/ricardo/projects/JavaProjects/src/ai/woman.jpg");
		BufferedImage image = null;
		
		image = ImageIO.read(file);
		
		System.out.println("Height  -> "+   image.getHeight());
		System.out.println("Weight -> "+ image.getWidth());
		
		KMeansImg tr = new KMeansImg(file);
	
		System.out.println("getPixels() ---> " + tr.getPixels().length);
		int [][] intVect = tr.getPixels();
		
		System.out.println(" intVect.length -> "+ intVect.length);
		System.out.println(" intVect[0].length -> " + intVect[0].length);
		
		/*for (int i=0; i < intVect.length; i++) {
			for (int j =0 ; j < intVect[0].length; j++) {
				System.out.println(intVect[i][j]);
				if (j == 2)
					break;
			}				
		}
		*/
		int xLenght = intVect.length;
		int yLenght = intVect[0].length;
		BufferedImage b = new BufferedImage(xLenght, yLenght, 3);

		for(int x = 0; x < xLenght; x++) {
		    for(int y = 0; y < yLenght; y++) {
		        int rgb = (int)intVect[x][y]<<16 | (int)intVect[x][y] << 8 | (int)intVect[x][y];
		        b.setRGB(y, x, rgb);
		    }
		}
		ImageIO.write(b, "jpg", new File("/home/ricardo/projects/JavaProjects/src/ai/Doublearray.jpg"));
		System.out.println("end");
		System.out.println(" image transformed in double \n");
		System.out.println( "img height -> " + b.getHeight());
		System.out.println(" img width -> "+ b.getWidth());
		
		int [][] imagx = new int[b.getHeight()][b.getWidth()];
		for (int i = 0; i < b.getHeight(); i++)
			for(int j=0; j < b.getWidth(); j++)
				imagx[i][j] = b.getRGB(i, j);

	/*System.out.println(" print imagx \n");	
		
	for (int i=0; i < imagx.length; i++) {
		for (int j =0 ; j < imagx[0].length; j++) {
			System.out.println(imagx[i][j]);
			if (j == 2)
				break;
		}				
	}
	*/

		Kmeans_var(imagx, 4, 100, 0.01);
		
	}

	public static void Kmeans_var(int [][] xin, int nc, int maxepoch, double tol) {
	
		double dJ = tol+1;
		int epoch =0;
		int [][] input = xin;
		double [] J = new double[maxepoch];
		int [] ndx = new int [xin.length];
		for(int i =0; i< ndx.length; i++) {
			ndx[i] = new Random().nextInt(nc);
		}
		
		int [][] c = new int[nc][xin[0].length];
		
		for (int n=0; n < nc; n++) {
			
			for (int t =0; t < c[0].length; t++) {
				
				c[n][t] = input[ndx[n]][t];
			//	System.out.println(c[n][t]);
			}
		}
		
		double [] dist = new double[nc];
		double sum =0;
		int um= 0;
		int [][] u = new int[xin.length][nc];
		
		while(dJ>tol && epoch <  maxepoch) {
			double [] dist1 = null ;
			for (int vector = 0; vector< xin.length; vector++) {
				
				for (int centro =0; centro < nc; centro++) {
							sum =0;
							
							for (int j =0; j< input[0].length; j++) {
								sum = sum + Math.pow((input[vector][j] - c[centro][j]), 2);
							}
							dist[centro] = sum;	
				}
				
				for (int i =0; i < dist.length; i++) {
					if ( dist[i] == min(dist)) {
						 um= i;
					}
				}
				
				u[vector][um] = 1;
				
			}
			int centro =0;
			for ( centro =0; centro < nc; centro++) {
				int [] tmp = new int [xin.length];
				for (int vector = 0; vector< xin.length; vector++) {
					tmp[vector] = u[vector][centro];
				}
				int [] unos = new int [xin.length];
				for (int i =0; i < tmp.length; i++) {
					if (tmp[i] == 1) {
						unos [i] =i;
						
					}
				}
				double soma1 =0;
				int [] soma2 = new int[xin[0].length];
				//continue
				for(int membro = 0; membro < unos.length; membro++) {
				//	soma1 = soma1 
					//soma1=soma1+sum(((input(unos(membro),:))-c(centro,:)).^2);
					double suma=0;
					soma1 = soma1+ suma;
					for (int j =0; j < input[0].length; j++) {
						
						suma = suma + Math.pow((input[unos[membro]][j] - c[centro][j]), 2);
						soma2[j] =  soma2[j] + (input[unos[membro]][j]);
					}
					
					 dist1 = new double[nc];
					
					dist1[centro] = soma1;
					
					for (int i =0; i < soma2.length; i++) {
						
						c[centro][i] = soma2[i] / unos.length;
					}
					
				}// membro
				
				
				
			}
			//epoch
		
		
			double sumd = 0;
			for (int i=0; i < dist1.length; i++) {
				sumd =  sumd + dist1[i];
			}
			J[epoch] = sumd;
			
			if (epoch ==1)
				dJ = J[epoch];
			else {
				if (epoch > 1)
					dJ = Math.abs(J[epoch-1] - J[epoch]);
			}
				
			System.out.println("epoch -> "+ epoch);
			
			epoch = epoch +1;
			
			
		}
		//-------
		
		int [][] saida = input;
		
		int [] clas= new int[xin.length];
		
		for (int i =0;  i< xin.length; i++) {
			int cla =0;
			for (int j=0; j< nc; j++) {
				
				if(u[i][j]!=0) {
					cla = i; 
				}
				
			}
			clas[i] = cla;
		}
		System.out.println("class " + clas.length);
	//	for (int i =0; i < clas.length; i++) {
		//	System.out.println(clas[i]);
			/*if (i ==30)
				break;*/
	//	}
		
	}// kmeans var
	
	public static double min(double[] vec) {
		/**
		 * Functions that finds the min val
		 * in a vector
		 * 
		 */
		double ele = vec[0];
		double min = 0;
		for (int x = 1; x < vec.length; x++) {
			double vectmp = vec[x];

			if (ele < vectmp) {
				min = ele;
			} else
				min = vectmp;

			ele = min;
		}

		return min;

	}

}
