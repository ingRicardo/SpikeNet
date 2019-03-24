/**
 * 
 */
package ai;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.util.ArrayList;

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
 * @author ricardo
 *
 */
public class TestKmeans extends JFrame{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final int NUM_CLUSTERS =2;
	private static final int TOTAL_DATA =7;
	private static final double [][] 	SAMPLES = {
						   {1.0, 1.0},
						   {1.5, 2.0},
						   {3.0, 4.0},
						   {5.0, 7.0},
						   {3.5, 5.0},
						   {4.5, 5.0},
						   {3.5, 4.5}};
	
	private static final String title = "Clustering k-means";
	private static ArrayList<Data> dataSet = new ArrayList<Data>();
	private static ArrayList<Centroid> centroids = new ArrayList<Centroid>();
	 XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
	 XYSeries NormalDataAdded = new XYSeries("NormalData");
	 XYSeries Class_A_Added = new XYSeries("Class_A");
	 XYSeries Class_B_Added = new XYSeries("Class_B");
	public TestKmeans(String s) {

		
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
	            	
	            	 xySeriesCollection.removeAllSeries();
	            	 
	            		// Print out clustering results.
	            	    for(int i = 0; i < NUM_CLUSTERS; i++)
	            	    {
	            	        System.out.println("Cluster " + i + " includes:");
	            	        for(int j = 0; j < TOTAL_DATA; j++)
	            	        {
	            	            if(dataSet.get(j).cluster() == i){
	            	                System.out.println("     (" + dataSet.get(j).X() + ", " + dataSet.get(j).Y() + ")");
	            	                if (i ==0) {
	            	                	Class_A_Added.add(dataSet.get(j).X(), dataSet.get(j).Y());
	            	                }
	            	                
	            	                if (i ==1) {
	            	                	Class_B_Added.add(dataSet.get(j).X(), dataSet.get(j).Y());
	            	                }
	            	                
	            	            }
	            	        } // j
	            	        System.out.println();
	            	    } // i
	            	 
	       	
	            	xySeriesCollection.addSeries(Class_A_Added);
	            	xySeriesCollection.addSeries(Class_B_Added);
	                
	            }
	        }));
	        
	        
	       
	        this.add(control, BorderLayout.SOUTH);
		
		
	}
	
private ChartPanel createNormalizedPanel() {
		
	initialize();
	kMeanCluster();
	
	// Print out clustering results.
    for(int i = 0; i < NUM_CLUSTERS; i++)
    {
        System.out.println("Cluster " + i + " includes:");
        for(int j = 0; j < TOTAL_DATA; j++)
        {
            if(dataSet.get(j).cluster() == i){
                System.out.println("     (" + dataSet.get(j).X() + ", " + dataSet.get(j).Y() + ")");
            }
        } // j
        System.out.println();
    } // i
    
    // Print out centroid results.
    System.out.println("Centroids finalized at:");
    for(int i = 0; i < NUM_CLUSTERS; i++)
    {
        System.out.println("     (" + centroids.get(i).X() + ", " + centroids.get(i).Y());
    }
    System.out.print("\n");
   

	
        JFreeChart jfreechart = ChartFactory.createScatterPlot(
            title, "X", "Y", createNormalizedData(),
            PlotOrientation.VERTICAL, true, true, false);
        XYPlot xyPlot = (XYPlot) jfreechart.getPlot();
        xyPlot.setDomainCrosshairVisible(true);
        xyPlot.setRangeCrosshairVisible(true);
        XYItemRenderer renderer = xyPlot.getRenderer();
        renderer.setSeriesPaint(0, Color.blue);
        NumberAxis domain = (NumberAxis) xyPlot.getDomainAxis();
        domain.setRange(1.00, 10.00);
        domain.setTickUnit(new NumberTickUnit(1.0));
        domain.setVerticalTickLabels(true);
        NumberAxis rangep = (NumberAxis) xyPlot.getRangeAxis();
        rangep.setRange(1.0, 10.0);
        rangep.setTickUnit(new NumberTickUnit(1.0));
        return new ChartPanel(jfreechart);
    
	}

private XYDataset createNormalizedData() {
    
    XYSeries NormalData = new XYSeries("NormalData");

    // Print out clustering results.
    for(int i = 0; i < NUM_CLUSTERS; i++)
    {
        System.out.println("Cluster " + i + " includes:");
        for(int j = 0; j < TOTAL_DATA; j++)
        {
            if(dataSet.get(j).cluster() == i){
                System.out.println("     (" + dataSet.get(j).X() + ", " + dataSet.get(j).Y() + ")");
             
                	NormalData.add(dataSet.get(j).X(), dataSet.get(j).Y());
               
                
                
                
            }
        } // j
        System.out.println();
    } // i
    
    xySeriesCollection.addSeries(NormalData);

    return xySeriesCollection;
}

	private static void initialize()
    {
        System.out.println("Centroids initialized at:");
        centroids.add(new Centroid(1.0, 1.0)); // lowest set.
        centroids.add(new Centroid(5.0, 7.0)); // highest set.
        System.out.println("     (" + centroids.get(0).X() + ", " + centroids.get(0).Y() + ")");
        System.out.println("     (" + centroids.get(1).X() + ", " + centroids.get(1).Y() + ")");
        System.out.print("\n");
        return;
    }
	
	
	public static void main(String[] args) {
		

		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	                TestKmeans demo = new TestKmeans(title);
	                demo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	                demo.pack();
	                demo.setLocationRelativeTo(null);
	                demo.setVisible(true);
	            }
	        });	
		
	}
	
	public static void test() {
		
		for (double[] ds : SAMPLES) {
			System.out.println(" ");
			printVec(ds);
		}
		
	}
	
	public static void printVec(double [] vec) {
	
		for(int i=0; i< vec.length; i++) {
			System.out.println(vec[i]);
		}
	}
	
	private static void kMeanCluster() {
		
		final double bigNumber = Math.pow(10, 10);
		double minimum = bigNumber;
		double distance = 0.0;
		int sampleNumber =0;
		int cluster =0;
		boolean isStillMoving = true;
		Data newData = null;
		
		while(dataSet.size() < TOTAL_DATA) {
			
			newData = new Data(SAMPLES[sampleNumber][0], SAMPLES[sampleNumber][1]);
			dataSet.add(newData);
			minimum = bigNumber;
			
			for (int i=0; i < NUM_CLUSTERS; i++) {
				distance = dist(newData, centroids.get(i));
				if(distance < minimum) {
					minimum =  distance;
					cluster =i;
				}
				
			}
			
			newData.cluster(cluster);
			
			for(int i=0; i < NUM_CLUSTERS; i++) {
				
				int totalX =0;
				int totalY =0;
				int totalInCluster =0;
				
				for (int j=0; j< dataSet.size(); j++) {
					
					if(dataSet.get(j).cluster() == i) {
						
						totalX += dataSet.get(j).X();
						totalY += dataSet.get(j).Y();
						totalInCluster++;
					}
					
				}
				if (totalInCluster >0) {
					centroids.get(i).X(totalX / totalInCluster);
					centroids.get(i).Y(totalY / totalInCluster);
					
				}
				
			}
			sampleNumber++;
			
			
		}
		
		while(isStillMoving) {
			
			for (int i=0; i< NUM_CLUSTERS; i++) {
				int totalX =0;
				int totalY =0;
				int totalInCluster =0;
				
				for (int j=0; j< dataSet.size(); j++) {
					
					if (dataSet.get(j).cluster() == i) {
						totalX += dataSet.get(j).X();
						totalY += dataSet.get(j).Y();
						totalInCluster++;
					}
					
				}
				if(totalInCluster > 0) {
					centroids.get(i).X(totalX / totalInCluster);
					centroids.get(i).Y(totalY / totalInCluster);
				}
				
			}
			
			isStillMoving = false;
			
			for (int i = 0; i < dataSet.size(); i++) {
				Data tempData = dataSet.get(i);
				minimum = bigNumber;
				
				for(int j=0; j < NUM_CLUSTERS; j++) {
					
					distance = dist(tempData, centroids.get(j));
					if (distance < minimum) {
						minimum = distance;
						cluster =j;
					}
				}
				
				tempData.cluster(cluster);
				if(tempData.cluster() != cluster) {
					tempData.cluster(cluster);
					isStillMoving = true;
				}
				
				
			}
			
		}
		return;
	}
	
	
	private static double dist(Data d, Centroid c) {
		return Math.sqrt(Math.pow(( c.Y() - d.Y()  ), 2) + Math.pow((c.X()-d.X()), 2));
	}
		
	private static class Data{
		
		private double mX = 0;
		private double mY = 0;
		private int mCluster =0;
		
		public Data(double x, double y)
		{
			this.X(x);
			this.Y(y);
			return;
		}
		public void X(double x) {
			
			this.mX = x;
			return;
			
		}
		public double X() {
			
			return this.mX;
		}
		
		public void Y(double y) {
			this.mY = y;
			return;
			
		}
		
		public double Y() {
			
			return this.mY;
		}
		
		public void cluster(int clusterNumber) {
			
			this.mCluster = clusterNumber;
			return;
		}
		
		public int cluster() {
			
			return this.mCluster;
		}
		
		
		
	}
	
	private static class Centroid{
		
		private double mX =0.0;
		private double mY = 0.0;
		
		public Centroid(double newX, double newY) {
			this.mX = newX;
			this.mY = newY;
			return;
		}
		
		public void X(double newX) {
			
			this.mX =  newX;
			return;
		}
		
		public double X() {
			
			return this.mX;
			
		}
		
		public void Y(double newY) {
			
			this.mY = newY;
			return;
		}
		public double Y() {
			
			return this.mY;
		}
			
		
	}
	
}
