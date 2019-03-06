/**
 * 
 */
package nnets;

/**
 * @author TIJUANA
 *
 */

import java.awt.*;       // Using AWT's Graphics and Color
import java.awt.event.*; // Using AWT event classes and listener interfaces
import java.util.ArrayList;

import javax.swing.*;    // Using Swing's components and containers
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


public class SpikeNetV3 extends JFrame{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private static final int N = 8;
    private static final String title = "Spike Net - Simple Clustering";
    //private static final Random rand = new Random();
    private XYSeries added = new XYSeries("Spikes");
	static double  ein[][] = {  {21.0,69.0,13.0,14.0,26.0,8.5,0.0,25.0,30.0},
								{3.0,7.0,19.0,4.0,15.0,5.0,9.8,1.0,2.0}};
	
	
	static double  ein_rf[][] = {  {21.0,69.0,13.0,14.0,26.0,8.5,0.0,25.0,30.0},
								{3.0,7.0,19.0,4.0,15.0,5.0,9.8,1.0,2.0}};
	
	
	static double[][] aus = ein;	
	static double[][] spike = new double[ein.length][ein[0].length];
	
	
	static double maxd = 0.2;

	boolean show = false;
	static int typ =1;
	
	static double [][] aus_rf;
	static double [][] cluster = null;
    public SpikeNetV3(String s) {
        super(s);
        final ChartPanel chartPanel = createDemoPanel();
        this.add(chartPanel, BorderLayout.CENTER);
        JPanel control = new JPanel();
        control.add(new JButton(new AbstractAction("Add Spikes") {

			private static final long serialVersionUID = 1L;

			@Override
            public void actionPerformed(ActionEvent e) {
            	 double movX =0.00;
            	
            	for (double[] sp : spike) {
                	double x = 0.0;
                    double y = 0.0; 
                    System.out.println(" ");
            		for (int i = 0; i < sp.length; i++) {
            			movX = movX + 0.035;
                        x = movX;
                        y = sp[i];                        
            			System.out.println("sp "+ y);            			
                        added.add(x, y);
                    }
				}            	                
            }
        }));
        this.add(control, BorderLayout.SOUTH);
    }

    private ChartPanel createDemoPanel() {
    	    	
    	kodieren_ln();
    	Kodieren_rf();
    	SimpleClustering();
    	
        JFreeChart jfreechart = ChartFactory.createScatterPlot(
            title, "X", "Y", createSampleData(),
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

    private XYDataset createSampleData() {
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        XYSeries series = new XYSeries("Delays");
       double movX =0.0;
        for (double[] ds : aus) {
        	double x = 0.0;
            double y = 0.0;        	
        	for (int i = 0; i < ds.length; i++) {
        		movX = movX + 0.035;
                x = movX;
                y = ds[i];
                series.add(x, y);
            }
		}
                
      //  xySeriesCollection.addSeries(series);
        xySeriesCollection.addSeries(added);
                
        XYSeries seriesRF = new XYSeries("Delays RF");
        movX =0;
        for (double[] dsrf : aus_rf) {
        	double x = 0.0;
            double y = 0.0;        	
        	for (int i = 0; i < dsrf.length; i++) {
        		movX = movX + 0.035;
                x = movX;
                y = dsrf[i];
                seriesRF.add(x, y);
            }
		}
     //   xySeriesCollection.addSeries(seriesRF);
        
        XYSeries series0 = new XYSeries("Class 0");
        XYSeries series1 = new XYSeries("Class 1");
        int r =0;
        movX =0;
        for (double[] ds : cluster) {
        	double x = 0.0;
            double y = 0.0;
        	for (int i =0; i< ds.length; i++) {
        		movX = movX + 0.035;
                x = movX;                
                y = ds[i];
                
        		if (r == 0)
        			series0.add(x,y);
        		
        		if (r == 1)
        			series1.add(x,y);
        		
        	}        	
			r++;        	        	        	
		}
        
        xySeriesCollection.addSeries(series0);
        xySeriesCollection.addSeries(series1);
        return xySeriesCollection;
    }

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	                SpikeNetV3 demo = new SpikeNetV3(title);
	                demo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	                demo.pack();
	                demo.setLocationRelativeTo(null);
	                demo.setVisible(true);
	            }
	        });
		
		
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
	
	public static void kodieren_ln() {
		
			double range =0.0;
			int r=0;
			
			System.out.println("print Data");
			for (double[] is : ein) {
				System.out.println(" ");
				print(is);
			}
			
			for (double[] row : ein) {
				
				range = MaxInArray(row) - MinInArray(row);
				double[] rown = new double[ein[0].length];
				for (int x =0; x < row.length; x++) {
					rown[x] = (row[x] - MinInArray(row))/range;	
				}
				ein[r] = rown;
				spike[r] = rown;
				r++;
			}
			
			System.out.println(" ");
			System.out.println("Normalized Data");
			for (double[] is : ein) {
				System.out.println(" ");
				print(is);
			}			
			 r=0;
			for (double[] row : aus) {
				
				double[] rown = new double[ein[0].length];
				
				for (int c=0; c < row.length; c++) {
					rown[c]= row[c]*maxd;
					
				}
			aus[r] = rown; 
				r++;	
			
			}
			System.out.println("");
			System.out.println("calculating delays");
			System.out.println("print AUS ");
					
			for (double[] ds : aus) {
				System.out.println("");				
				print(ds);
				
			}	    			
	}
	
	public static void Kodieren_rf() {
		
		System.out.println(" Kodieren_rf");
		
		int[] rf =new int[ein_rf.length];
		int i =0;
		
		double rho = 0.5;
		double[] sigma = new double [rf.length];
						
		double f_cut = 0.5;
		
		for (double[] ds : ein_rf) {
			System.out.println("");
			rf[i] = 1;
			print(ds);
			i++;
			
		}
		System.out.println("RF ->> "+ rf);
		
		for (i=0; i < rf.length; i++) {
			System.out.println( rf[i]);
		}
		
		int t_neur = rf.length;		
		System.out.println(" t_neur ->> "+ t_neur);
		int t_cj = ein[0].length;
		System.out.println(" t_cj ->> "+ t_cj);
		int t_in = rf.length;
		System.out.println(" t_in ->> "+ t_in);
		double [][] aus = new double[ t_cj][t_in];	
		double [][] zeros = new double[ t_cj][t_in];
		
		System.out.println(" \nFILL MATRIX WITH ZEROS");
		int con =0;
		
		for (double[] is : zeros) {	
			System.out.println(" ");
			for (int x =0 ; x < is.length; x++) {				
				zeros [con][x] =	0.0;
			}
			print(zeros[con]);
			con++;
			
		}
		
		aus = zeros;
		
		System.out.println(" \n FILL AUS MATRIX WITH ZEROS");
		con =0;
		
		for (double[] is : aus) {	
			System.out.println(" ");
			for (int x =0 ; x < is.length; x++) {				
				aus [con][x] =	0.0;
			}
			print(aus[con]);
			con++;
			
		}
		
		double[] ct = new double[t_in]; 
		
		for (i =0 ; i < t_in ; i++ )
			ct[i] = 0.0; 
		
		System.out.println("PRINT CT ");
		
		print(ct);
				
		if (rf.length < 3) {
			System.out.println(" Receptive Field must be greater than 3");
	
		}
		
		double range;
		int r =0;
		for (double[] row : ein) {
			
			range = MaxInArray(row) - MinInArray(row);
			double[] rown = new double[ein_rf[0].length];
			for (int x =0; x < row.length; x++) {
				rown[x] = (row[x] - MinInArray(row))/range;	
			}
			ein_rf[r] = rown;
			r++;
		}
		r=0;
		
		System.out.println("PRINT ein_rf normalized");
		
		for (double[] d : ein_rf) {
			System.out.println(" ");
			print(d);
		}
		
		double  cdis = 0.0;
		int j =0;
		if (typ == 1) {
			for (j =1; j < rf.length; j++) {
				cdis = 1/rf[j];
				ct[j] = (j-0.5)*cdis;
				sigma [j] = 1 / (1.5*(t_in -2));
			}
		}else if (typ ==0){
			for (j =1; j < rf.length; j++) {
				cdis = 1/(rf[j]-2);
				ct[j] = (j-1.5)*cdis;
				sigma [j] = 1 / (1.5*(t_in +1));
			}
						
		}else {
			System.out.println("Invalid typ parameter it needs to be 1 or 0");			
		}
		
		System.out.println("\n PRINT distance between centers ");
		
		print(ct);
		
		System.out.println("\n PRINT SIGMAS ");
		print(sigma);
		
		double aux =0.0;
		int n=0;
		
		System.out.println("spike threshold ->> " + maxd*f_cut);
		for (j=1 ; j < t_cj; j++) {
			System.out.println("");
			for (int k =1 ; k < rf.length; k++) {
				aux = maxd - maxd * Math.exp( (-1*(( Math.pow((ein_rf[k][j]- ct[k]), 2) ) / (2*( Math.pow(sigma[k], 2) )))));
			
				aus[j][k] = (aux * (1/rho))* rho;
												
				if (aus[j][k] > maxd*f_cut) {
					aus[j][k] = -1;
				}
				System.out.println("aus ->> " + aus[j+n][k]);
			}
			
		}
		
		aus_rf = aus;
		
	}
	
	public static void SimpleClustering() {
	  
	  System.out.println("\nSimpleClustering");
		
	  double  aux1[][] =ein;		 
      double rho = 0.1, tau =1.84, teta=0.2;
	  int maxd = 1, tmax=8;		
	  int tam=0, in_neu =0 , out_neu = 2;

	  for (double[] ds : aux1) {
		in_neu++;
		tam =0;
		for (int x =0; x < ds.length; x++)
			tam++;
	  }
	  in_neu = in_neu +1 ;
	  System.out.println("\ntam "+ tam + "\nin_neu "+ in_neu + "\nout_neu "+ out_neu);				
	  int tra = tam; double trb = (tam/2)+ tra;
	  System.out.println("trb ->> "+ trb);		
	  System.out.println("* PRINT AUX1 *");
		
	  double [][] aus = kodieren_ln(aux1, maxd, rho, 0);						
	  double [][] ca  = { {2.0, 2.5},
						  {7.8, 6.5}};
	  double sum =0.0;
	  int c=0;
	  for (double[] ds : aus) {		
		  sum = 0.0;
		  for (int i=0; i<tra ; i++) {
			  sum = sum + ds[i];
		  }						   
		  ca[c][0] = ((sum/tra)/rho)*rho;
		  ca[c][1] = ((sum/tra)/rho)*rho + 0.05;		 
		  c++;
	  }
	  	  	 
	  double [][] cb  = { {8.2, 6.4},
			  			  {3, 3.2} };				
	  c=0;
	  for (double[] ds : aus) {		
		  sum = 0.0;
		  for (int i=0; i<tra ; i++) {
			  sum = sum + ds[i];
		  }						   
		  cb[c][0] = ((sum/trb)/rho)*rho;
		  cb[c][1] = ((sum/trb)/rho)*rho + 0.18;
		 
		  c++;
	  }
		  
	  double maxdm = 3;	  
	  double ca1 = ca[0][0] , ca2 = ca[0][1] ;
	  double cb1 = cb[0][0] , cb2 = cb[0][1] ;	  	  
	  double [][] delaytmp = {{ca1, cb1},
			  				  {ca2, cb2},
			  				  { maxdm,  maxdm }};
	  
	  System.out.println( "\nca1 : "+ ca1 + " cb1 : "+ cb1 );
	  System.out.println( "ca2 : "+ ca2 + " cb2 : "+ cb2 );
	  
	  double [][] d = new double[delaytmp.length][delaytmp[0].length];
		  
		int le=0; 
		  for (double[] ds : delaytmp) {
			  
			  for (int x =0; x < ds.length; x++) {
				  d[le][x] = maxdm - ds[x];				  				
			  }
		  
		   le++;
		}
		  
		  System.out.println("\nPRINT DELAYS d");
		  for (double[] es : d) {
			  System.out.println(" ");
			  print(es);
			
		}
		  
		  double w[][] = {	{1,1.42},
				  			{1,1.38},
				  			{1,1}};
		  		  
		  
		  double [][] ausa = new double [aus.length + 1][aus[0].length];		  
		  System.out.println("\nTesting ausa = aus \n");
		  
		  int r=0;
		  for (double[] es : ausa) {
			
			  if (r < aus.length)
				  for(int x = 0; x < es.length; x++) {
					  ausa[r][x] = aus[r][x];
				  
				  }
			  else {
				  for(int x = 0; x < es.length; x++) {
					  ausa[r][x] = 1 * maxdm;
				  
				  }
			  }
			  r++;			 				  
		}
		  		  

		  System.out.println("AUGMENTED MATRIX AUS");
		  
		  for (double[] e : ausa) {
			  System.out.println(" ");
			print(e);
		}
		 	 
		 double  neu ;
		 double t; 
		 double [] out = null;		  
		 ArrayList<Double> ltmp = new ArrayList<Double>();		
		
		 double dt =0.0, sai=0.0;
		 for (int k=0 ; k <tam ; k++) {
			 t=0;
			 neu = 0;
			 while (neu == 0 && t < tmax ) {			
				out =  new double[out_neu]; 				
					for (int j =0; j < out_neu; j++) {			
						for (int i =0; i < in_neu; i++) {							
							dt = t- ausa[i][k];									
							sai = w[i][j] * (((dt - d[i][j]) /tau) * (Math.exp(1 - ((dt - d[i][j])/tau) )));
							if (sai < 0)
								sai = 0;
							
							out[j] = out[j] + sai;															
						}
						
					}										
					if (MaxInArray(out) >= teta) {														
							for (int x = 0; x < out.length; x++) {
								 if (out[x] == MaxInArray(out)) {
									 neu = x;
									 ltmp.add(neu);									 									 
								 }																	
							}													 													
					}
			
					t = t +rho;
			 }			

		 }
		  		 	
		 System.out.println("\nltmp -- size "+ ltmp.size());
		 
		 int tot = (in_neu-1)*  aus[0].length;
		 System.out.println("tot ->> "+ tot);
		 
		 int mo =0;
		 if (tot != ltmp.size()) {
			 mo = tot%ltmp.size();
			 System.out.println("mod ->> "+  mo );
		 }
		 int i =0;
		 for (i =0 ; i < mo; i++) {			 
			 ltmp.add(0.0);
		 }
		
		 System.out.println("\nQuantity of O ADDED : "+ i+ " \n");
		 
		double [][] auscla = new double [aus.length][aus[0].length];
		
		r=0;
		int s =0;
		for (double[] es : auscla) {
			
			for (i=0; i< es.length; i++) {
				auscla[r][i] = ltmp.get(s);
				s++;
			}
			r++;	
			
		}		 		
		 cluster = new double [out_neu][tam];
	     int ro=0;
		 
		 ArrayList<Double> clas0 = new ArrayList<Double>();
		 ArrayList<Double> clas1 = new ArrayList<Double>();
		 		 		 
		 for (i =0; i< out_neu; i++) {
			 
			 for (double[] es : auscla) {
				
				for (int j =0; j< es.length; j++) {
					 for (i =0; i< out_neu; i++) {
						if (i == auscla[ro][j] ) {
							
							if (auscla[ro][j] == 1) {
								System.out.println("class "+auscla[ro][j] + " val -> "+ aus[ro][j]);						
								clas1.add(aus[ro][j]);
							}
														
							if (auscla[ro][j] == 0) {
								System.out.println("class "+auscla[ro][j] + " val -> "+ aus[ro][j]);								
								clas0.add(aus[ro][j]);
								
							}
													
						 }
					}
						
				}
				 ro++;
			}
			
		 }
		 
		 r =0;
		 for (double[] row : cluster) {
			
			 for (int x =0; x < row.length; x++) {
				 if (r == 0)
					 cluster[r][x] = clas0.get(x);
				 if(r == 1)
					 cluster[r][x] = clas1.get(x);
				 
			 }
			 r++;
			 
		}
		 
		 System.out.println("\n<-- PRINT CLUSTERS -->");
		 c =0;
		 for (double[] es : cluster) {
			 System.out.println("\nclass -> "+c);
			 for (i =0; i< es.length; i++) {
				 System.out.println(es[i]);
			 }
			c++; 
		}
		  
	}

	
	
	public static double[][] kodieren_ln(double [][] aux1, int maxd, double rho, int plt){

		double range =0.0;
		int r=0;
		
		System.out.println("print Data");
		for (double[] is : aux1) {
			System.out.println(" ");
			print(is);
		}
		
		for (double[] row : aux1) {
			
			range = MaxInArray(row) - MinInArray(row);
			double[] rown = new double[aux1[0].length];
			for (int x =0; x < row.length; x++) {
				rown[x] = (row[x] - MinInArray(row))/range;	
			}
			aux1[r] = rown;
			r++;
		}
		
		System.out.println(" ");
		System.out.println("Normalized Data");
		for (double[] is : aux1) {
			System.out.println(" ");
			print(is);
		}
		
		 r=0;
		for (double[] row : aux1) {
			
			double[] rown = new double[aux1[0].length];
			
			for (int c=0; c < row.length; c++) {
				rown[c]= row[c]*maxd;			
			}
			aux1[r] = rown; 
			r++;	
		
		}

		return aux1;
	}
	
}
