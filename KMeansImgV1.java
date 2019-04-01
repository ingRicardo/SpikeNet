/**
 * 
 */
package nnets;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.Random;

import javax.imageio.ImageIO;

/**
 * @author TIJUANA
 *
 */
public class KMeansImgV1 {

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
	
	static int [] clas = null;
	static int [][] c = null;
	
	public KMeansImgV1(File file) {
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
		File file = new File("C:/projects/JavaProjects/src/nnets/woman.jpg");
		BufferedImage image = null;
		
		image = ImageIO.read(file);
		
		
		
		System.out.println("Height  -> "+   image.getHeight());
		System.out.println("Weight -> "+ image.getWidth());
		
		KMeansImgV1 tr = new KMeansImgV1(file);
	
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
		ImageIO.write(b, "jpg", new File("C:/projects/JavaProjects/src/nnets/Doublearray.jpg"));
		System.out.println("end");
		System.out.println(" image transformed in double \n");
		System.out.println( "img height -> " + b.getHeight());
		System.out.println(" img width -> "+ b.getWidth());
		
		int [][] imagx = new int[b.getHeight()][b.getWidth()];
		for (int i = 0; i < b.getHeight(); i++)
			for(int j=0; j < b.getWidth(); j++)
				imagx[i][j] = b.getRGB(i, j);
		
	//	int [][] imagx = imageToPixels(image);
		
	/*	for (int i =0; i < imagx.length; i++) {
			for (int j =0; j < imagx[0].length; j++) {
				System.out.println(" imagx[ "+i+" ][ "+j+" ] -> " + imagx[i][j]);
			}
		}
		*/
		
	/*	 xLenght = imagx.length;
		 yLenght = imagx[0].length;
		 b = new BufferedImage(xLenght, yLenght, 3);

		for(int x = 0; x < xLenght; x++) {
		    for(int y = 0; y < yLenght; y++) {
		        int rgb = (int)imagx[x][y]<<16 | (int)imagx[x][y] << 8 | (int)imagx[x][y];
		        b.setRGB(x, y, rgb);
		    }
		}
		ImageIO.write(b, "jpg", new File("C:/projects/JavaProjects/src/nnets/imagx.jpg"));
		System.out.println("end");
		
		*/
	/*System.out.println(" print imagx \n");	
		
	for (int i=0; i < imagx.length; i++) {
		for (int j =0 ; j < imagx[0].length; j++) {
			System.out.println(imagx[i][j]);
			if (j == 2)
				break;
		}				
	}
	*/
		int nc = 4;
		Kmeans_var(imagx, nc, 100, 0.01);
		//int [][] imagx = new int[b.getHeight()][b.getWidth()]; 
		// imagx =  xin
		//nc = 4;
		//int [][] c = new int[nc][xin[0].length];
		//int [] clas= new int[xin.length];
		// imagy[i] = clas[i] = b.getHeight()
	//	clas= new int[b.getHeight()];
		int [][]  c_int = c;
		int [][] imagy = new int[imagx.length][c_int[0].length];
		
		for (int i=0; i < imagx.length; i++) {
			
				for (int j =0; j < imagy[0].length; j++) {
					for (int x =0; x < c_int.length; x++) {
						if (clas[i] == x) {
							imagy[i][j]	= c_int[x][j];
						}
					}
				}
		}
		
		
		System.out.println(" <- imagy -> ");
		
		/*System.out.println(imagy.length);
		System.out.println(imagy[0].length);
		for (int  i=0; i < c_int.length; i++) {
			System.out.println("");
			for (int j =0; j < c_int[0].length; j++) {
				System.out.println("imagy[ "+i+" ][ "+j+" ] -> "+ imagy[i][j] );
			
			}
		}*/
		/////////////////////////////
		
		xLenght = imagy.length;
		yLenght = imagy[0].length;
		BufferedImage imagz = new BufferedImage(xLenght, yLenght, 3);

		for(int x = 0; x < xLenght; x++) {
		    for(int y = 0; y < yLenght; y++) {
		        int rgb = imagy[x][y]<<16 | imagy[x][y] << 8 | imagy[x][y];
		        imagz.setRGB(y, x, rgb);
		    }
		}
		ImageIO.write(imagz, "jpg", new File("C:/projects/JavaProjects/src/nnets/Doublearray2.jpg"));
		System.out.println("end");
		
		
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
		
	/*	for (int i=0; i < ndx.length; i++ ) {
			
			System.out.println("ndx "+ ndx[i]);
		}
		System.out.println();
		*/
		c = new int[nc][xin[0].length];
		
		for (int i=0; i < input.length; i++)
			for (int n=0; n < nc; n++) {
				
				
				for (int t =0; t < c[0].length; t++) {
					
				//	if (ndx[n] == input[i][t] ) {
						c[n][t] = input[ndx[n]][t];
				//	}

				//	System.out.println(c[n][t]);
				}
			}
		
	/*	System.out.println("< -c ->");
		System.out.println("row -> "+ c.length);
		System.out.println("col -> "+ c[0].length);
		for (int i=0; i < c.length; i++) {
			System.out.println("");
			for (int j =0; j < c[i].length; j++)
				System.out.println(c[i][j]);
		}
		
		*/
		double [] dist = new double[nc];
		double sum =0;
		int um= 0;
		int [][] u = new int[xin.length][nc];
		
		while(dJ>tol && epoch <  maxepoch) {
			double [] dist1 = new double[nc];
			for (int vector = 0; vector< xin.length; vector++) {
				
				for (int centro =0; centro < nc; centro++) {
							sum =0;
							
							for (int j =0; j< input[0].length; j++) {
								sum = sum + Math.pow((input[vector][j] - c[centro][j]), 2);
							}
							dist[centro] = sum;	
				}
				
				for (int i =0; i < dist.length; i++) {
				//	System.out.println("dist -> "+ dist[i]);
					if ( dist[i] == min(dist)) {
						
						 um= i;
					//	 System.out.println("um -> " + um);
					}
				}
				
				u[vector][um] = 1;
				
			}
			int centro =0;
			for ( centro =0; centro < nc; centro++) {
				
				/*
				 * for centro=1:nc % para cada centro:
					unos=find(u(:,centro)==1); % acha vetores da classe nc
					soma1=0; % zera somas
					soma2=zeros(1,col);
					for membro=1:size(unos,1) % soma distancias da classe nc
					soma1=soma1+sum(((input(unos(membro),:))-c(centro,:)).^2);
					soma2=soma2+(input(unos(membro),:));
					end
					105
				 * 
				 */
				int [] unos = new int [xin.length];
				for (int vector = 0; vector< xin.length; vector++) {
					if (u[vector][centro] ==1) {
						unos [vector] =vector;
					//	System.out.println(" unos[ "+vector+" ] -> " +   unos[vector]);
					}
				}
				
		//		System.out.println(" ");

				double soma1 =0;
				int [] soma2 = new int[xin[0].length];
				//continue
				for(int membro = 0; membro < unos.length; membro++) {
				//	soma1 = soma1 
					//soma1=soma1+sum(((input(unos(membro),:))-c(centro,:)).^2);
					double suma=0;
					
					for (int j =0; j < input[0].length; j++) {
						
						suma = suma + Math.pow((input[unos[membro]][j] - c[centro][j]), 2);
						soma2[j] =  soma2[j] + (input[unos[membro]][j]);
						c[centro][j] = soma2[j] / unos.length;
					}
					soma1 = soma1+ suma;
				//	System.out.println(" soma1 -> " + soma1);
					dist1[centro] = soma1;
					
		//		System.out.println("dist1[ "+centro+" ] --> "+ dist1[centro]);	
					
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
		
		clas= new int[xin.length];
		
		int cla =0;
		
		for (int i =0;  i< xin.length; i++) {
		
			for (int j=0; j< nc; j++) {
				
				if(u[i][j]!=0) {
					cla = i; 
				}
				
			}
			clas[i] = cla;
		}
		/*
		System.out.println(" u mat ");
		System.out.println(u.length);
		System.out.println(u[0].length);
		
		for (int i=0; i < u.length; i++) {
			System.out.println(" ");
			for(int j =0; j < u[0].length; j++) {
				System.out.println(u[i][j]);
			}
		}
		*/
		/*	System.out.println("class " + clas.length);
		for (int i =0; i < clas.length; i++)
			System.out.println(clas[i]);*/
			
		
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