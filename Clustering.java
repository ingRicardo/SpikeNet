/**
 * 
 */
package nnets;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * @author I.C. Ricardo Alfredo Macias Olvera
 *
 */

public class Clustering extends JFrame{
	/**
	 * Clustering with SNN (Spike Neural Network)
	 * Multiple synapses architecture
	 * 5 clusters with gaussian distribution with means 0
	 * Each cluster has 3 points (data)
	 * output classes : 5
	 * codification : fixed receptive fields
	 * Number of receptive fields : 1
	 * Type of synapses : multiple
	 */
	private static final long serialVersionUID = 1L;
	private static final String title = "Clustering SNN";
//	private static int out_neu=5;
	 XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
	 XYSeries NormalDataAdded = new XYSeries("NormalData");
	 XYSeries DelaysAdded = new XYSeries("Delays");
	 
	 XYSeries Class_A_Added = new XYSeries("Class_A");
	 XYSeries Class_B_Added = new XYSeries("Class_B");
	 XYSeries Class_C_Added = new XYSeries("Class_C");
	 XYSeries Class_D_Added = new XYSeries("Class_D");
	 XYSeries Class_E_Added = new XYSeries("Class_E");
		// Create a new ArrayList 
    static ArrayList<Integer> newListClas = new ArrayList<Integer>(); 
	static  ArrayList<Integer> arClas= new ArrayList<Integer>();

	static double [][] aus = null;
	private static  double [][] ein = null;
	
	static double [] vclasA = null;
	static double [] vclasB = null;
	static double [] vclasC = null;
	static double [] vclasD = null;
	static double [] vclasE = null;
	
	
	public Clustering(String s) {

		
		final ChartPanel chart_Normal_Panel = createNormalizedPanel();
		this.add(chart_Normal_Panel, BorderLayout.CENTER);
		
		 JPanel control = new JPanel();
	        control.add(new JButton(new AbstractAction("Clear") {
	        	/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

				@Override
	            public void actionPerformed(ActionEvent e) {
	            	
	        		xySeriesCollection.removeAllSeries();
	            	            	                
	            }
	        }));
	        
	        
	        
	        control.add(new JButton(new AbstractAction("Add Clusters") {

				private static final long serialVersionUID = 1L;

				@Override
	            public void actionPerformed(ActionEvent e) {
	            	 double movX =0.00;
	            	
	            	for (double sp : vclasA) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp;                        
	            			System.out.println("A -> "+ y);            			
	            			Class_A_Added.add(x, y);
	                   
					} 
	            	
	            	
	            	
	            	for (double sp : vclasB) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp;                        
	            			System.out.println("B ->"+ y);            			
	            			Class_B_Added.add(x, y);
	                    
					}
	            	
	            	
	            	for (double sp : vclasC) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp;                        
	            			System.out.println("C ->"+ y);            			
	            			Class_C_Added.add(x, y);
	                    
					}
	            	
	            	
	            	for (double sp : vclasD) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp;                        
	            			System.out.println("D ->"+ y);            			
	            			Class_D_Added.add(x, y);
	                    
					}
	            	
	            	for (double sp : vclasE) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp;                        
	            			System.out.println("E ->"+ y);            			
	            			Class_E_Added.add(x, y);
	                   
					}
	            	
	            	xySeriesCollection.addSeries(Class_E_Added);
	            	xySeriesCollection.addSeries(Class_D_Added);
	            	xySeriesCollection.addSeries(Class_C_Added);
	            	xySeriesCollection.addSeries(Class_B_Added);
	            	xySeriesCollection.addSeries(Class_A_Added);          	                
	            }
	        }));
	        
	        
	        
	        control.add(new JButton(new AbstractAction("Add Normal") {

				private static final long serialVersionUID = 1L;

				@Override
	            public void actionPerformed(ActionEvent e) {
	            	 double movX =0.00;
	            	
	            	for (double[] sp : ein) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		for (int i = 0; i < sp.length; i++) {
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp[i];                        
	            			System.out.println("sp "+ y);            			
	            			NormalDataAdded.add(x, y);
	                    }
					} 
	            	  xySeriesCollection.addSeries(NormalDataAdded);          	                
	            }
	        }));
	        
	        
	        
	        
	        
	        
	        control.add(new JButton(new AbstractAction("Add Delays") {

				private static final long serialVersionUID = 1L;

				@Override
	            public void actionPerformed(ActionEvent e) {
	            	 double movX =0.00;
	            	
	            	for (double[] sp : aus) {
	                	double x = 0.0;
	                    double y = 0.0; 
	                    System.out.println(" ");
	            		for (int i = 0; i < sp.length; i++) {
	            			movX = movX + 0.035;
	                        x = movX;
	                        y = sp[i];                        
	            			System.out.println("d "+ y);            			
	            			DelaysAdded.add(x, y);
	                    }
					} 
	            	  xySeriesCollection.addSeries(DelaysAdded);          	                
	            }
	        }));
	        
	        
	        
	        
	        this.add(control, BorderLayout.SOUTH);
		
		
	}
	
	private ChartPanel createNormalizedPanel() {
		
		test();		
        JFreeChart jfreechart = ChartFactory.createScatterPlot(
            title, "X", "Y", createNormalizedData(),
            PlotOrientation.VERTICAL, true, true, false);
        XYPlot xyPlot = (XYPlot) jfreechart.getPlot();
        xyPlot.setDomainCrosshairVisible(true);
        xyPlot.setRangeCrosshairVisible(true);
        XYItemRenderer renderer = xyPlot.getRenderer();
        renderer.setSeriesPaint(0, Color.blue);
        NumberAxis domain = (NumberAxis) xyPlot.getDomainAxis();
        domain.setRange(0.00, 1.00);
        domain.setTickUnit(new NumberTickUnit(0.1));
        domain.setVerticalTickLabels(true);
        NumberAxis rangep = (NumberAxis) xyPlot.getRangeAxis();
        rangep.setRange(0.0, 1.0);
        rangep.setTickUnit(new NumberTickUnit(0.1));
        return new ChartPanel(jfreechart);
    
	}
	
	
    private XYDataset createNormalizedData() {
       
        double movX =0.0;
        XYSeries NormalData = new XYSeries("NormalData");
    
        movX =0;
        for (double[] ds : ein) {
        	double x = 0.0;
            double y = 0.0;
        	for (int i =0; i< ds.length; i++) {
        		movX = movX + 0.035;
                x = movX;                
                y = ds[i];
                    
                NormalData.add(x,y);        	        	        	
        	}
        
        }
        
        xySeriesCollection.addSeries(NormalData);
    
        return xySeriesCollection;
    }
				
	/**
	 * @param args
	 */
	public static void main(String[] args) {

		// TODO Auto-generated method stub
		
		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	                Clustering demo = new Clustering(title);
	                demo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	                demo.pack();
	                demo.setLocationRelativeTo(null);
	                demo.setVisible(true);
	            }
	        });								
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
		
		/***
		 * 
		 *  Encode analog values by means of multiple fixed receptive fields 
		 */
		
		int t_neur = 0;
		for (int i =0; i < rf.length; i++) {
			t_neur = t_neur + rf[i]; 
			
		}
		
		int t_cj =  ein.length;
		
		System.out.println( "\n t_cj - > t_cj -> t_cj -> "+ t_cj);
		
		aus = new double [t_cj][ t_neur]; 
				
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
				
			}
		}
		
					
		return aus;
	}
	public static void test() {
		
		int out_neu=5;
		double rho = 0.6;
		int conj = 5;
		int tau = 4;
		int [] rf = {1, 1, 1};
		double [] sigma = {1/(1.5*rf[0]-2),1/(1.5*rf[1]-2), 1/(1.5*rf[0]-2)};
		double f_cut =0.9;
		int maxd = 2;
		
	
		int [] d = {0,1,2};
		int w_max = 1, w_min=0;
		int max_epoch = 2;
		int t_max =5;
		double eta = 0.25;
		double beta = 0.1;
		double nu = 5.0;
	
		double kapa = 1 - (Math.pow(nu, 2)) / (2* Math.log(beta/(beta+1)) ); 
		
		int in_neu = 0;
		for (int i =0 ; i < rf.length; i++) {
			in_neu = in_neu + rf[i];
			
		}
		
		System.out.println("in_neu -> "+ in_neu);
		int ssin = d.length;
			
		//double teta = 1.5* (in_neu/out_neu) * (w_max-w_min)*ssin/2;
		double teta = 0.5;
		// Initialize weights
			
		double [][][] w= new double[in_neu][out_neu][ssin];
		
		double [][] ein = {{2.3,6, 8.5,4.5},
						   {9.3,3.3, 4.2,7.3},
						   {15,3.5,7.2,8.2},
						   {22.1,6.6, 5.7, 4.4},
						   {1.3,8.4,1.6,1}};
		Clustering.ein = ein;
		
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
		aus = kodieren_rf(ein, rf, maxd, sigma, f_cut, rho, 0, 0);
		
		System.out.println("\n aus after kodieren_rf ------->  "+ aus.length);
		
		int ctrl_1=0;
		System.out.println("\n delays \n");
		
		for (double [] d1 : aus) {
			System.out.println(" ");
			print(d1);
		}
		
		double [][][] delta_t = new double [1][1][ssin];
		
		System.out.println(" <<<<<<<<  delta_t - ssin >>>>>>>>>> "+ ssin);
		
		double [][] f_times = new double[2][1];
		
		
		System.out.println(" \n\n Training begins\n");
		int neu = 0;
		double inst =0;
		double dneu =0;
		
		while(ctrl_1 <= max_epoch) {
			
			System.out.println(" ");
			
			for (double[] tr : aus) {
				
				System.out.println(" ");
					for (int z = 0 ; z < tr.length; z ++) {
					
						System.out.print (" "+ tr[z] );
						
						ArrayList<ArrayList<Object>> firings = wer_schiesst(tr, t_max, w, d, tau, rho, teta, 0);
						
						for(int i=0; i< firings.get(0).size(); i++) {								
							neu =  (int) firings.get(0).get(i);						
							dneu = neu;								
							f_times[0][i] =  dneu;																			
						}
						for (int j =0 ; j < firings.get(1).size(); j++) {							
							f_times[1][j] = (double) firings.get(1).get(j);
							inst = (double) firings.get(1).get(j);
						}
						
						if ( neu!=0 && firings.get(0).size() == 1  ) { //size of neu							
							for (int i=0; i< in_neu; i++) {
																					
								 if (tr[z] == -1) {

									 for ( i=0; i < w[z][neu].length; i++) {
										 w[z][neu][i] = w[z][neu][i] -eta;

									 }
								 }else {
									 System.out.println("< -delta_t -> ");									 
									 for ( i=0; i < ssin; i++) {
										
										 delta_t[0][0][i] = tr[z] +d[z] -inst;										 
 w[z][neu][i] =  w[z][neu][i] + eta * ((1 +beta) * Math.exp(-1*( (Math.pow((delta_t[0][0][i]) +3.3, 2))/(2*kapa-2)))-beta);
 

									 }

									 
								 }
								 
								 System.out.println("\nGetting indx < w_min  \n");
								 for ( i=0; i < w[z][neu].length; i++) {
									 if (w[z][neu][i] < w_min) {
										 System.out.println(" w_min  i -> -> ->  "+ i);
										 w[z][neu][i] = w_min;
									 }else if (w[z][neu][i] > w_max) {
										 System.out.println(" w_max i -> -> ->  "+ i);
										 w[z][neu][i] = w_max;
									 }
									 									 										
								 }
								 
							}
						}
					}
						
			}
			
			ctrl_1++;
			teta = teta + (0.3*teta) /max_epoch;
		}
		int i=0, j=0, ausr =0, ausc=0;
				
		for ( i=0; i < aus.length; i++) {
			ausc = 0;
			for (j =0; j < aus[i].length; j++) {
				ausc = ausc +1 ;
			}
			ausr = ausr + 1;
		}
		System.out.println("aus groups -> "+ ausr);
		System.out.println("aus data -> "+ ausc);		
		System.out.println(" conj -> "+ conj);
		System.out.println(" in_neu -> "+ in_neu);
		
		
		//Testing the net
		int [] clas = new int[conj];
		
		int c =0;
		for (c=0; c < conj; c++) {
			for (double[] tr : aus) {
					
					ArrayList<ArrayList<Object>> firings = wer_schiesst(tr, t_max, w, d, tau, rho, teta, 0);
					
					for(j=0; j< firings.get(0).size(); j++) {								
						neu =  (int) firings.get(0).get(j);	
						System.out.println("neu ---testing--->  " + neu);															
					}

			}
			
			clas[c] = neu;
		}
		
		System.out.println("\n PRINT CLASES -> \n ");
		for(int x =0; x< clas.length; x++)
			System.out.println(clas[x]);
		
		System.out.println(" ");
		
		System.out.println("total clases arr - > " + arClas.size());
		   
        for (Integer element : arClas) {   
            if (!newListClas.contains(element)) {   
            	newListClas.add(element); 
            } 
        } 
  
		System.out.println("NO DUPLICATES CLASES");
		
		for (Integer integer : newListClas) {
			System.out.println(" clas -----> "+ integer);
		}
		
		
		 vclasA = new double[ein[0].length];
		 vclasB = new double[ein[0].length];
		 vclasC = new double[ein[0].length];
		 vclasD = new double[ein[0].length];
		 vclasE = new double[ein[0].length];
		
		for ( i=0; i< out_neu; i++) {
			
			for ( j=0; j< newListClas.size(); j++) {
				if (i == newListClas.get(j)) {
					
					System.out.println("ein[ "+ i +" ] -> "+ ein[i]);
					
					for (int z =0; z < ein[i].length; z++) {
						
						System.out.println("val -> " + ein[i][z]);
						
						if (i == newListClas.get(j) && i == 0) {							
							vclasA[z] = ein[i][z];
						}						
						if (i == newListClas.get(j) && i == 1) {							
							vclasB[z] = ein[i][z];
						}						
						if (i == newListClas.get(j) && i == 2) {							
							vclasC[z] = ein[i][z];
						}						
						if (i == newListClas.get(j) && i == 3) {							
							vclasD[z] = ein[i][z];
						}						
						if (i == newListClas.get(j) && i == 4) {							
							vclasE[z] = ein[i][z];
						}						
					}																
				}
			}
			
		}
		
		System.out.println(" size of vectors -> \n");
		
		System.out.println("\n vclasA.length -> "+ vclasA.length);		
		print(vclasA);
		System.out.println("\n vclasB.length -> "+ vclasB.length);		
		print(vclasB);
		System.out.println("\n vclasC.length -> "+ vclasC.length);		
		print(vclasC);
		System.out.println("\n vclasD.length -> "+ vclasD.length);
		print(vclasD);
		System.out.println("\n vclasE.length -> "+ vclasE.length);
		print(vclasE);
		
	}
	
	
	
	public static ArrayList<ArrayList<Object>> wer_schiesst(double[] in,int t_max, double [][][] w, int [] d, int tau, double rho, double teta, int plt ) {
		
		/**
		 * Function that returns the neuron was fired and the firing time 
		 * 
		 */
				
		System.out.println(" ");
		System.out.println("in -> "+ in+  " t_max -> " + t_max +" w -> "+ w + " d -> "+ d + " tau -> "+ tau+ " "
							+ " rho -> "+ rho + " teta -> "+ teta + " plt -> "+ plt );
		
		int ein = 0, ausx =0, ssin=0;
		for (double [][] input : w) {
			ausx=0;
			for (double[] output : input) {
				ssin=0;
				for (double ssyn : output) {					
					ssin++;
				}
				ausx++;				
			}			
		ein++;
		}
		System.out.println("Wiegths dimensions");
		System.out.println("\ninput - > "+ ein + " output -> "+ ausx + " ssyn -> "+ ssin);
		// input - > 3 output -> 5 ssyn -> 3
				
		double [] indx = new double[in.length];
		for (int i=0; i < in.length; i++) {
			
			if(in[i] !=-1) {
				indx[i] = in[i];
			}
		}
						
		int tam = indx.length;
		System.out.println("indx  ->> "+ tam);
		System.out.println("\n <- indx -> \n");
		print(indx);
		
		int neu =0;
		double inst =0;
		
		ArrayList<ArrayList<Object>> firings = new ArrayList<ArrayList<Object>>();
		
		double [] vals = new double[2];
		if (indx.length == 0) {
			neu =0; inst = -1;
			vals[0] = neu;
			vals[1] = inst;
			return firings;
		}
		
		System.out.println("\n d size "+ d.length + " ssin "+ ssin);
	
		 
		double [] out ;
		
		double dt=0, dtt=0, sai=0,t =0;
		
		ArrayList<Object> neuIndx = null;
		ArrayList<Object> neuInst = null;
		
		System.out.println("j -> "+ aus[0].length);
		System.out.println(" i -> "+ tam);
		System.out.println(" k -> "+ ssin);
		
		//ArrayList<Integer> indxTemp = new ArrayList<Integer>();
		while (neu == 0 && t <= t_max) {
			
			System.out.println("<- wer_schiesst -> ");
			
			out = new double[ausx];
			
			for (int j=0; j< ausx; j++) {
								
				System.out.println(" ");
				for (int i=0; i < tam; i++) {
					
					dt = t - indx[i];											
					for (int k =0; k < ssin; k++) {
						
						dtt= (dt - d[k])/tau;																								
						sai = w[i][j][k]*dtt*Math.exp(1-dtt);
												
						if (sai < 0)
							sai =0;

						out[j] = out[j] +sai;

					}

				}
				
			}
			
			if (max(out)>= teta ) {
				neuInst = new ArrayList<Object>();
				neuIndx = new ArrayList<Object>();

				for (int x =0; x < out.length ; x++) {

					if (out[x] == max(out)) {
						neu = x;
						neuIndx.add(neu);	
						arClas.add(neu);
						inst = t;
						System.out.println("neu - > "+ neu + " inst -> "+ inst);
						neuInst.add(inst);
												
					}else
						neu = 0;
											
				}
				firings.add(neuIndx);
				firings.add(neuInst);					
			}else
				neu =0;
			t=t+rho;	
		}
					
		return firings;
	}

}
