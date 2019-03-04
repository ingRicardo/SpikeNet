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


public class SpikeNetV2 extends JFrame{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//private static final int N = 8;
    private static final String title = "Spike Net";
    //private static final Random rand = new Random();
    private XYSeries added = new XYSeries("Spikes");
	static double  ein[][] = {{21.0,19.0,13.0,14.0,26.0,28.0,0.0,95.0,20.0},
			   				  {8.0,7.0,14.0,4.0,15.0,5.0,3.0,1.0,2.0},
			   				  {43.6, 5.5, 1.0, 33.4, 87.0, 22.0, 17.0, 44.0, 8.9}};
	
	
	static double  ein_rf[][] = {{21.0,19.0,13.0,14.0,26.0,28.0,0.0,95.0,20.0},
			   					 {8.0,7.0,14.0,4.0,15.0,5.0,3.0,1.0,2.0},
			   					{43.6, 5.5, 1.0, 33.4, 87.0, 22.0, 17.0, 44.0, 8.9}};
	
	
	static double[][] aus = ein;	
	static double[][] spike = new double[ein.length][ein[0].length];
	
	
	static double maxd = 0.2;

	boolean show = false;
	static int typ =1;
	
	static double [][] aus_rf;
	
    public SpikeNetV2(String s) {
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
        
        
        xySeriesCollection.addSeries(series);
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
        xySeriesCollection.addSeries(seriesRF);
        
        return xySeriesCollection;
    }

    
    
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	                SpikeNetV2 demo = new SpikeNetV2(title);
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
					//System.out.println(rown[c]);
					
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

}
