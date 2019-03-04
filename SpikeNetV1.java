/**
 * 
 */
package nnets;

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

/**
 * @author RICA
 *
 */
public class SpikeNetV1 extends JFrame {

	private static final long serialVersionUID = 1L;
	//private static final int N = 8;
    private static final String title = "Spike Net";
    //private static final Random rand = new Random();
    private XYSeries added = new XYSeries("Spikes");
	static double  ein[][] = {{21.0,19.0,13.0,14.0,26.0,28.0,0.0,95.0,20.0},
			   {8.0,7.0,14.0,4.0,15.0,5.0,3.0,1.0,2.0}};
	static double[][] aus = ein;	
	static double[][] spike = new double[ein.length][ein[0].length];
    
    public SpikeNetV1(String s) {
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
            			movX = movX + 0.05;
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
    	
    	
		double maxd = 0.2;
	//	double rho = 0.5;
		//boolean show = false;
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
        		movX = movX + 0.05;
                x = movX;
                y = ds[i];
                series.add(x, y);
            }
		}
        
        
        xySeriesCollection.addSeries(series);
        xySeriesCollection.addSeries(added);
        return xySeriesCollection;
    }

	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	                SpikeNetV1 demo = new SpikeNetV1(title);
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
	
	public static void kodieren_ln(int ein[][], double maxd, double rho, boolean plot) {
		
		
		
	}
	
}




