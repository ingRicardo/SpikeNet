/**
 * 
 */
package nnets;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.EventQueue;
import java.util.ArrayList;

import javax.swing.JFrame;

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
 * @author I.C Ricardo Alfredo Macias Olvera - TIJUANA
 * 	 *  - Clustering using SNN
	 *  -It uses 2 groups with Gaussian distribution with zero as mean
	 *    Each group has 8 Random - values
	 *  -It uses a method to learn in a supervised way with
	 *   axonal delays. The weights are equals to 1 
	 *   
	 *   -Out classes : 2 
	 *   -Codification : Linear with one neuron as a reference
	 *   -Receptive Fields : 0
	 *   -Type of Synapses : simples
 */
public class SimpleClusterV2 extends JFrame{
	private static final long serialVersionUID = 1L;

	/**
	 * @param args
	 */

	private static final int NUM_CLUSTERS = 2;

	private static final double[][] SAMPLES = { 
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
												{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 } };
	
	private static final int TOTAL_DATA = SAMPLES.length; //Total Data
	
	static double[][] N_SAMPLES = null; 
	private static ArrayList<Data> dataSet = new ArrayList<Data>();
	static ArrayList<DataA> dataSetAumented = null;
	static int tam; // number of points
	static int in_neu; // number of input neurons
	static int out_neu; // number of output neurons
	static double rho; // time
	static double tau; // Time EPSP
	static int maxd; // encoded interval
	static int tmax; // max interval
	static double teta;
 
	private static final String title = "Simple Cluster SNN";
	XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
	//XYSeries Delays = new XYSeries("DelaysData");
	
	XYSeries Class_A_Added = new XYSeries("Class_A");
	XYSeries Class_B_Added = new XYSeries("Class_B");
	 
	public SimpleClusterV2(String title) {
		final ChartPanel chart_Delay_Panel = createDelayPanel();
		this.add(chart_Delay_Panel, BorderLayout.CENTER);
	}
	private ChartPanel createDelayPanel() {
		
	
        JFreeChart jfreechart = ChartFactory.createScatterPlot(
            title, "X", "Y", createDelayData(),
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
	
    private XYDataset createDelayData() {
            
     // Print out clustering results.
	    for(int i = 0; i < NUM_CLUSTERS; i++)
	    {

	        for(int j = 0; j < TOTAL_DATA; j++)
	        {
	            if(dataSetAumented.get(j).cluster() == i){

	                if (i ==0) {
	                	Class_A_Added.add(dataSetAumented.get(j).X(), dataSetAumented.get(j).Y());
	                }
	                
	                if (i ==1) {
	                	Class_B_Added.add(dataSetAumented.get(j).X(), dataSetAumented.get(j).Y());
	                }
	                
	            }
	        } // j
	        System.out.println();
	    } // i
	 
	    xySeriesCollection.addSeries(Class_A_Added);
	    xySeriesCollection.addSeries(Class_B_Added);
   
    
        return xySeriesCollection;
    }
	
	public static void Initialize() {
		tam = SAMPLES.length;
		in_neu = SAMPLES[0].length + 1;
		out_neu = 2;
		rho = 0.1;
		tau = 1.84;
		maxd = 10;
		tmax = 6;
		teta = 1.5;

	}

	public static void Training() {
		int tra = 2;
		int trb = (tam / 2) + tra;
		double ca, cb;
		ArrayList<Data> aus = kodieren_ln(); // encoding inputs
		
	
		
		double sum = 0;
		for (int i = 0; i < tra; i++) {
			sum = sum + (aus.get(i).X() + aus.get(i).Y());

		}
		System.out.println(" sum -> "+ sum);
		ca = Math.round(((sum / tra) / rho) * rho); // means of class 1
		System.out.println("\n ca -> " + ca);
		System.out.println("trb -> " + trb);

		sum = 0;
		// change here
		for (int i = tra +1; i < trb; i++) {
			sum = sum + (aus.get(i).X() + aus.get(i).Y());
		}
		cb = Math.round(((sum / tra) / rho) * rho); // means of class 2
		System.out.println("\n cb -> " + cb);

		double[] d = { maxd - ca, maxd - cb, maxd - maxd };

		for (int i = 0; i < d.length; i++) {
			System.out.println(" d[ " + i + " ] -> " + d[i]);

		}

		double[] w = { 1, 1, 1 };
		int sampleNumber = 0;

		dataSetAumented = new ArrayList<DataA>();
		DataA newData = null;
		System.out.println(" aus.size() -> " + aus.size());
		while (dataSetAumented.size() < TOTAL_DATA) {
			newData = new DataA(N_SAMPLES[sampleNumber][0], N_SAMPLES[sampleNumber][1], 1 *maxd);
			dataSetAumented.add(newData);
			sampleNumber++;
		}
		System.out.println("\n dataSetAumented " + dataSetAumented.size() + " tam->  " + tam + "\n");
		for (int j = 0; j < dataSetAumented.size(); j++)
			System.out.println("X -> " + dataSetAumented.get(j).X() + " Y -> " + dataSetAumented.get(j).Y() + " Z -> "
					+ dataSetAumented.get(j).Z());

		double dt = 0, sai = 0;
		double t = 0;
		double[] out = null; // new double[out_neu];	
		int neu=0;	
		for (int k = 0; k < dataSetAumented.size(); k++) {
			t = 0;
			neu=0;
			while ( neu == 0 &&  t <= tmax) {
				out = new double[in_neu];
				for (int j = 0; j < NUM_CLUSTERS; j++) {
					for (int i = 0; i < in_neu; i++) {

						if (i == 0)
							dt = t - dataSetAumented.get(k).X();

						if (i == 1)
							dt = t - dataSetAumented.get(k).Y();

						if (i == 2)
							dt = t - dataSetAumented.get(k).Z();

						sai = w[i] * ((dt - d[i]) / tau) * Math.exp(1 - ((dt - d[i]) / tau));

						if (sai < 0)
							sai = 0;
											
						out[i] = out[i] + sai;
						System.out.println("out[ " + i + " ] - > " + out[i]);
					} // in_neu
					System.out.println("");

				} // out neu
				
				double max = max(out);
				System.out.println(" max -> "+ max);
				if (max >= teta) {				
					for (int to = 0; to < out.length ; to++) {
						if (out[to] == max) {							 
							neu = to;
							
						}
					}
					
				} // max

				t = t + rho;
			} // while
			System.out.println("neu -> "+ neu);
			dataSetAumented.get(k).cluster(neu);
		} // tam
		
		sampleNumber = 0;

		for (int i = 0; i < NUM_CLUSTERS; i++) {

			System.out.println("Cluster " + i + " includes: ");

			for (int j = 0; j < TOTAL_DATA; j++) {

				if (dataSetAumented.get(j).cluster() == i) {
					System.out.println(" (" + dataSetAumented.get(j).X() + ", " + dataSetAumented.get(j).Y() + ")");
				}
			}
			System.out.println();

		}

	}

	// linear encoding
	public static ArrayList<Data> kodieren_ln() {

		int sampleNumber = 0;
		Data newData = null;
		Data normalData = null;

		while (dataSet.size() < TOTAL_DATA) {
			newData = new Data(SAMPLES[sampleNumber][0], SAMPLES[sampleNumber][1]);
			dataSet.add(newData);
			sampleNumber++;
		}

		double[] vecX = new double[dataSet.size()];
		double[] vecY = new double[dataSet.size()];

		double[] range = new double[SAMPLES[0].length];
		double[] max = new double[SAMPLES[0].length];
		double[] min = new double[SAMPLES[0].length];

		N_SAMPLES = new double[SAMPLES.length][SAMPLES[0].length];

		for (int s = 0; s < SAMPLES[0].length; s++)
			for (int i = 0; i < dataSet.size(); i++) {

				if (s == 0)
					vecX[i] = dataSet.get(i).X();

				if (s == 1)
					vecY[i] = dataSet.get(i).Y();

			}

		for (int s = 0; s < SAMPLES[0].length; s++) {

			if (s == 0) {
				max[s] = max(vecX);
				min[s] = min(vecX);
				range[s] = max[s] - min[s];
			}

			if (s == 1) {
				max[s] = max(vecY);
				min[s] = min(vecY);
				range[s] = max[s] - min[s];
			}

		}

		for (int s = 0; s < SAMPLES[0].length; s++) {
			for (int i = 0; i < dataSet.size(); i++) {
				if (s == 0)
					N_SAMPLES[i][s] = (vecX[i] - min[s]) / range[s];
				if (s == 1)
					N_SAMPLES[i][s] = (vecY[i] - min[s]) / range[s];
			}
		}

		sampleNumber = 0;

		ArrayList<Data> dataSetNormal = new ArrayList<Data>();
		while (dataSetNormal.size() < TOTAL_DATA) {
			normalData = new Data(N_SAMPLES[sampleNumber][0], N_SAMPLES[sampleNumber][1]);
			dataSetNormal.add(normalData);
			sampleNumber++;
		}

		System.out.println("\nNormal dataSet [0-1]\n");
		for (int j = 0; j < dataSetNormal.size(); j++)
			System.out.println("X -> " + dataSetNormal.get(j).X() + " Y -> " + dataSetNormal.get(j).Y());

		ArrayList<Data> dataSetDelays = dataSetNormal;
		for (int j = 0; j < dataSetNormal.size(); j++) {
			dataSetDelays.get(j).X(Math.round(((dataSetNormal.get(j).X() * maxd) / rho) * rho));
			dataSetDelays.get(j).Y(Math.round(((dataSetNormal.get(j).Y() * maxd) / rho) * rho));
		}

		System.out.println("\n dataSetDelays \n");
		for (int j = 0; j < dataSetDelays.size(); j++)
			System.out.println("X -> " + dataSetDelays.get(j).X() + " Y -> " + dataSetDelays.get(j).Y());

		return dataSetDelays;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		Initialize();
		Training();
		
		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	            	SimpleClusterV2 demo = new SimpleClusterV2(title);
	                demo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	                demo.pack();
	                demo.setLocationRelativeTo(null);
	                demo.setVisible(true);
	            }
	        });		

	}

	private static class Data {

		private double mX = 0;
		private double mY = 0;
		private int mCluster = 0;

		public Data(double x, double y) {
			this.X(x);
			this.Y(y);
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

	private static class DataA {

		private double mX = 0;
		private double mY = 0;
		private double mZ = 0;
		private int mCluster = 0;

		public DataA(double x, double y, double z) {
			this.X(x);
			this.Y(y);
			this.Z(z);
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

		public void Z(double z) {
			this.mZ = z;
			return;
		}

		public double Z() {
			return this.mZ;
		}

		public void cluster(int clusterNumber) {

			this.mCluster = clusterNumber;
			return;
		}

		public int cluster() {
			return this.mCluster;
		}

	}

	public static double max(double[] vec) {
		double ele = vec[0];
		double max = 0;
		for (int x = 1; x < vec.length; x++) {
			double vectmp = vec[x];

			if (ele > vectmp) {
				max = ele;
			} else
				max = vectmp;

			ele = max;
		}

		return max;
	}

	public static double min(double[] vec) {

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
