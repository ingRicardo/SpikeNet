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
 * @author I.C. Ricardo Alfredo Macias Olvera TIJUANA
 *	/**
	 * Clustering with SNN (Spike Neural Network)
	 * SPIKE RESPONSE MODEL (SRM)  
	 * Multiple synapses architecture
	 * 5 clusters with gaussian distribution with means 0
	 * Sample has 10 data
	 * output classes : 5
	 * codification : fixed receptive fields
	 * Number of receptive fields : 1
	 * Type of synapses : multiple
	 */
 
public class ClusteringV2x extends JFrame{
	private static final long serialVersionUID = 1L;
	private static ArrayList<Data> dataSet = new ArrayList<Data>();
	private static ArrayList<Weight> weightSet = new ArrayList<Weight>();
	private static ArrayList<Data> NormalDataSet = new ArrayList<Data>();
	private static ArrayList<Data> DelayDataSet = new ArrayList<Data>();	
	private static final int NUM_CLUSTERS = 5;
	
	private static final double[][] SAMPLES = { 
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((10.0-1.0)+1))+ 1.0, (Math.random()*((10.0-1.0)+1))+ 1.0 } };
	
	private static final int TOTAL_DATA = SAMPLES.length;
	private static double rho; //time
	private static int tau; // EPSP time
	private static int [] rf; // number of neurons for input
	private static double [] sigma; //field amplitude
	private static double f_cut;
	private static int maxd; //encoded interval
	
	//--------------------------------------------
	private static int [] d;// sub synapses delays
	private static int w_max, w_min; // max and min weight
	private static int max_epoch; // max number of epochs
	private static int t_max; //max time of training
	private static double eta; // learning rate
	private static double beta; // learning window
	private static double nu;// neighborhood
	
	//----------------
	private static double kapa; // number for learning window
	private static int ssin ; // number of sub-synapses
	private static double teta;// threshold
	
	//----------------
	private static int typ; // curve center
	private static double inst;
	private static int neu;
	
	
	private static final String title = "Clustering SNN";	
	XYSeriesCollection xySeriesCollection = new XYSeriesCollection();		
	XYSeries Class_A_Added = new XYSeries("Class_A");
	XYSeries Class_B_Added = new XYSeries("Class_B");
	XYSeries Class_C_Added = new XYSeries("Class_C");
	XYSeries Class_D_Added = new XYSeries("Class_D");
	XYSeries Class_E_Added = new XYSeries("Class_E");
	
	public ClusteringV2x(String title) {
		final ChartPanel chart_Delay_Panel = createPanel();
		this.add(chart_Delay_Panel, BorderLayout.CENTER);
	}
	private ChartPanel createPanel() {
		
		
        JFreeChart jfreechart = ChartFactory.createScatterPlot(
            title, "X", "Y", createData(),
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
	
	
    private XYDataset createData() {
        
        // Print out clustering results.
   	    for(int i = 0; i < NUM_CLUSTERS; i++)
   	    {
   	        for(int j = 0; j < TOTAL_DATA; j++)
   	        {
   	            if(NormalDataSet.get(j).cluster() == i){

   	                if (i ==0) {
   	                	Class_A_Added.add(NormalDataSet.get(j).X(), NormalDataSet.get(j).Y());
   	                }
   	                
   	                if (i ==1) {
   	                	Class_B_Added.add(NormalDataSet.get(j).X(), NormalDataSet.get(j).Y());
   	                }
   	                
   	                if (i ==2) {
   	                	Class_C_Added.add(NormalDataSet.get(j).X(), NormalDataSet.get(j).Y());
   	                }
   	                
   	                if (i ==3) {
   	                	Class_D_Added.add(NormalDataSet.get(j).X(), NormalDataSet.get(j).Y());
   	                }
   	                
   	                if (i ==4) {
   	                	Class_E_Added.add(NormalDataSet.get(j).X(), NormalDataSet.get(j).Y());
   	                }
   	                
   	            }
   	        } // j
   	        System.out.println();
   	    } // i
   	 
   	    xySeriesCollection.addSeries(Class_A_Added);
   	    xySeriesCollection.addSeries(Class_B_Added);
   	    xySeriesCollection.addSeries(Class_C_Added);
   	    xySeriesCollection.addSeries(Class_D_Added);
   	    xySeriesCollection.addSeries(Class_E_Added);
   	    
           return xySeriesCollection;
       }
	
	
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		initialize();
		Training();
		
		 EventQueue.invokeLater(new Runnable() {

	            @Override
	            public void run() {
	            	ClusteringV2x demo = new ClusteringV2x(title);
	                demo.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	                demo.pack();
	                demo.setLocationRelativeTo(null);
	                demo.setVisible(true);
	            }
	        });		
	}
	
	public static void initialize() {
		rho = 0.1;
		tau =8;
	    rf = new int[2];
		rf[0] = 3;
		rf[1] = 3;
		sigma = new double[2];
		sigma[0] = 1/(1.5*(rf[0]-2));
		sigma[1] = 1/(1.5*(rf[1]-2));
		f_cut=0.3;
		maxd =5;
		//---------------------------
		d = new int[] { 0, 1, 2, 3};
		w_max=1; w_min=0;
		max_epoch=3;
		t_max = 10;
		eta = 0.35;
		beta = 0.2;
		nu =5.0;
		//-------------------------
		kapa = 1 - Math.pow(nu, 2) / (2 * Math.log(beta / (beta+1)));
		ssin = d.length;
		teta = 6; //(1.5*(in_neu/out_neu)*(w_max-w_min)*ssin/2); //12
		//----------------------------
		typ =0;
		
		int sampleNumber =0;
		
		Data newData = null;
		while (dataSet.size() < TOTAL_DATA) {
			newData = new Data(SAMPLES[sampleNumber][0], SAMPLES[sampleNumber][1]);
			dataSet.add(newData);
			sampleNumber ++;
		}
		
		//-----------------------------
		int wieghtSample =0;
		Weight newWeight = null;
		Ssin newSsin = null;
		while(weightSet.size() < TOTAL_DATA) {
			newSsin = new Ssin(w_min + ((Math.random()*((1.0-0)+1))+ 0)* (w_max -w_min),
					w_min + ((Math.random()*((1.0-0)+1))+ 0)* (w_max -w_min),
					w_min + ((Math.random()*((1.0-0)+1))+ 0)* (w_max -w_min),
					w_min + ((Math.random()*((1.0-0)+1))+ 0)* (w_max -w_min));
			
			newWeight = new Weight(
					w_min + ((Math.random()*((1.0-0)+1))+ 0)* (w_max -w_min),					
					w_min + ((Math.random()*((1.0-0)+1))+ 0)* (w_max -w_min),
					newSsin);
			
			weightSet.add(newWeight);
			
			wieghtSample++;
		}
		

	}
	
	public static void Training() {
		
		/**
		 * Training the network 
		 * updating weights , delta_t and teta
		 */
				
	kodieren_rf(); // encode input [0-1]
		
	int ctrl_1 =0;
	while ( ctrl_1 <= max_epoch) { // epochs
		
		double [] delta_t = new double[ssin];  // time motion
		for (int i =0; i< DelayDataSet.size(); i++) {
			wer_schiesst_(i); // firing 
		if (neu !=0) {
							
				if (DelayDataSet.get(i).X() ==-1 || DelayDataSet.get(i).Y() == -1) {
					
					if (neu == i) {
						
						weightSet.get(neu).X(weightSet.get(neu).X()-eta);
						weightSet.get(neu).Y(weightSet.get(neu).X()-eta);
						
						weightSet.get(neu).ssin.X(  weightSet.get(neu).ssin.X() -eta );
						weightSet.get(neu).ssin.Y(  weightSet.get(neu).ssin.Y() -eta );
						weightSet.get(neu).ssin.Z(  weightSet.get(neu).ssin.Z() -eta );
						weightSet.get(neu).ssin.W(  weightSet.get(neu).ssin.W() -eta );
					}
						
				}else {
					
					if (neu == i) {
						// firing neuron index
						
						for (int p =0; p < ssin; p++) { /* delta_t updating*/
							delta_t[p] = (DelayDataSet.get(i).X() +  DelayDataSet.get(i).Y()) + d[p] -inst;
						}
						// weights updating
						weightSet.get(neu).X(
											weightSet.get(neu).X()+ eta * 
												( (1 +beta) *
												  Math.exp(-1 * (Math.pow((delta_t[0]+2.3), 2) / (2*kapa-2) ))-beta));
						
						weightSet.get(neu).Y(
								weightSet.get(neu).Y()+ eta * 
									( (1 +beta) *
									  Math.exp(-1 * (Math.pow((delta_t[0]+2.3), 2) / (2*kapa-2) ))-beta));
						
						// sub- synapses updating
						weightSet.get(neu).ssin.X(
								weightSet.get(neu).ssin.X() + eta * 
									( (1 +beta) *
									  Math.exp(-1 * (Math.pow((delta_t[0]+2.3), 2) / (2*kapa-2) ))-beta));
						
						weightSet.get(neu).ssin.Y(
								weightSet.get(neu).ssin.Y() + eta * 
									( (1 +beta) *
									  Math.exp(-1 * (Math.pow((delta_t[1]+2.3), 2) / (2*kapa-2) ))-beta));
						weightSet.get(neu).ssin.Z(
								weightSet.get(neu).ssin.Z() + eta * 
									( (1 +beta) *
									  Math.exp(-1 * (Math.pow((delta_t[2]+2.3), 2) / (2*kapa-2) ))-beta));
						weightSet.get(neu).ssin.W(
								weightSet.get(neu).ssin.W() + eta * 
									( (1 +beta) *
									  Math.exp(-1 * (Math.pow((delta_t[3]+2.3), 2) / (2*kapa-2) ))-beta));
												
					
					
						//weight low limit must be minor or equals than w_min
						
						if (weightSet.get(neu).ssin.X() < w_min) {
							weightSet.get(neu).ssin.X(w_min);						
						}
						if (weightSet.get(neu).ssin.Y() < w_min) {
							weightSet.get(neu).ssin.Y(w_min);
						}
						if (weightSet.get(neu).ssin.Z() < w_min) {
							weightSet.get(neu).ssin.Z(w_min);
						}
						if (weightSet.get(neu).ssin.W() < w_min) {
							weightSet.get(neu).ssin.W(w_min);
						}
						
						//weight high limit must be greater or equals than w_max
						
						if (weightSet.get(neu).ssin.X() > w_max) {
							weightSet.get(neu).ssin.X( w_max);						
						}
						if (weightSet.get(neu).ssin.Y() > w_max) {
							weightSet.get(neu).ssin.Y(w_max);
						}
						if (weightSet.get(neu).ssin.Z() > w_max) {
							weightSet.get(neu).ssin.Z(w_max);
						}
						if (weightSet.get(neu).ssin.W() > w_max) {
							weightSet.get(neu).ssin.W(w_max);
						}
					
					}//neu
					
				
					
				}// else
				
			} // !=0
			
		}// for
		
		ctrl_1++; // epoch counting
		
		teta = teta+ (0.3*teta) / max_epoch; // updating 
	 }//while
	
	//-------------test the net------------------
	
	System.out.println(" <------ test the net ----->");
	for (int i=0; i< NormalDataSet.size(); i++) {
		
		wer_schiesst_( i );
				
		if (neu == i) {
			NormalDataSet.get(i).cluster(neu);
		}
		
	}
	

	for (int i = 0; i < NUM_CLUSTERS; i++) {

		System.out.println("Cluster " + i + " includes: ");

		for (int j = 0; j < TOTAL_DATA; j++) {

			if (NormalDataSet.get(j).cluster() == i) {
				System.out.println(" (" + NormalDataSet.get(j).X() + ", " + NormalDataSet.get(j).Y() + ")");
			}
		}
		System.out.println();

	}
	
}
	public static void wer_schiesst_(int ix ) {
		

		int [] indX = new int[DelayDataSet.size()];
		int [] indY = new int[DelayDataSet.size()];
		for (int z =0; z< DelayDataSet.size(); z++) {
			
			if ( DelayDataSet.get(z).X() != -1 ) { 

				indX[z] = z; // getting indexes of fired neurons
			}else {
				indX[z] = 100;
			}

			if ( DelayDataSet.get(z).Y() != -1 ) {

				indY[z] = z; // getting indexes of fired neurons
			}else {
				indY[z] = 100;
			}
	
		}
		
		double t =0, dt=0,dtt = 0, sai=0;			
		neu =0;
			while (neu == 0 && t <= t_max ) { // max time
			
			double [] out = new double[DelayDataSet.size()];
							
				if (indX[ix] == ix) { // data index
					
					dt = t - (DelayDataSet.get(ix).X() + DelayDataSet.get(ix).Y()); // time
					
					for (int k = 0; k < ssin; k++ ) { // sub- synapses
																
						dtt = (dt - d[k]) / tau; // delta time
																		
	/*SRM */				sai = 	weightSet.get(ix).X() * dtt * Math.exp(1 - dtt) +
								weightSet.get(ix).Y() * dtt * Math.exp(1 - dtt) +						
								weightSet.get(ix).ssin.X() * dtt * Math.exp(1 - dtt) +
								weightSet.get(ix).ssin.Y() * dtt * Math.exp(1 - dtt) +
								weightSet.get(ix).ssin.Z() * dtt * Math.exp(1 - dtt) +
								weightSet.get(ix).ssin.W()* dtt * Math.exp(1 - dtt);
				
						if (sai < 0)
							sai =0;
						
						out[ix] = out[ix] +sai;	

					}//ssin
						
				}//indx

			if (max(out) >= teta) { // threshold
				for (int j=0; j< out.length; j++) {
					
					if (max(out) == out[j]) {
						neu = j; // firing neuron index
						inst = t; // firing time
					}	
				}
				
			}
			t = t + rho; // time moving
			}			
	}
	

	public static void kodieren_rf() {
	/**
	 *  - function to normalized data [0-1]
	 *  - calculate delays
	 */
		if ( rf[0] >= 3 && rf[1] >= 3) { // input neurons size
			
			int [] ct = {rf[0], rf[1]}; //centers 
			
			double [] tmpVec = new double[SAMPLES[0].length];
			double [] tmpVecX = new double[SAMPLES[0].length];
			double [] tmpVecY = new double[SAMPLES[0].length];
			double range, min, max;
		
			int j =0;
			for ( j =0; j< SAMPLES[0].length; j++) {
				tmpVec = new double[SAMPLES.length];
				for (int i =0; i < SAMPLES.length; i++) {
					tmpVec[i] =SAMPLES[i][j];
				}
				min = min(tmpVec);
				max = max(tmpVec);
				range= max - min;
				
				for (int i =0; i < tmpVec.length; i++) {
					tmpVec[i]=(tmpVec[i] - min)/range;
				}
							
				if ( j == 0)
					tmpVecX = new double[SAMPLES.length];
				
				if (j == 1)
					tmpVecY = new double[SAMPLES.length];
				
				for (int i =0; i < tmpVec.length; i++) {
					
					if ( j == 0)
						tmpVecX[i] = 	tmpVec[i];
					
					if ( j == 1)
						tmpVecY[i] = 	tmpVec[i];
										
				}
				
			}
			Data newData = null;
			int numberSample =0;
			while(NormalDataSet.size() < TOTAL_DATA) { // data normalized [0-1]
				
				newData = new Data(tmpVecX[numberSample], tmpVecY[numberSample]);
				NormalDataSet.add(newData);
				numberSample++;
			}
			

		double cdis =0;  // distance between centers
		if (typ == 1) {
			
			cdis = 1/rf[0];
			
		} else if ( typ == 0) {
			
			cdis = 1/(rf[0]-2);
		}
		
		for (int i =0; i< rf.length; i++) { 
			
			ct[i] = (int) ((i -0.5) * cdis); // centers values
		}
		
		// getting delays 
		double auxX =0, auxY =0;
		double [] ausX = new double[TOTAL_DATA] ;
		double [] ausY = new double[TOTAL_DATA] ;

		Data newDelay = null;
			
		for (int i =0; i < NormalDataSet.size(); i++) { // computing delays
	
					 auxX = maxd - maxd *
							 Math.exp(-1 * (Math.pow((NormalDataSet.get(i).X() -ct[0]), 2)
									 / (2* Math.pow(sigma[0], 2))));
 
					 auxY = maxd - maxd *
							 Math.exp(-1 * (Math.pow((NormalDataSet.get(i).Y() -ct[1]), 2)
									 / (2* Math.pow(sigma[1], 2))));

				
					 ausX[i] = Math.round(auxX*(1/rho)) * rho; // round
					 ausY[i] = Math.round(auxY*(1/rho)) * rho; // round
				 
				 if ((ausX[i] > ( maxd * f_cut)) && (ausY[i] > ( maxd * f_cut))) { // firing or not

					 ausX[i] = -1;			// delay of -1
					 ausY[i] = -1;
				 }
				 
			
			newDelay = new Data(ausX[i], ausY[i]);		
			DelayDataSet.add(newDelay);
	
			
		}
		 
			
		}else {
			System.out.println("Incorrect size of receptive fields it must be greater than 2");
		}
		
	}
	
	
	private static class Ssin { // sub synapses class

		private double mX = 0;
		private double mY = 0;
		private double mZ = 0;
		private double mW = 0;

		public Ssin(double x, double y) {
			this.X(x);
			this.Y(y);
		}

		public Ssin(double x, double y, double z, double w) {
			this.X(x);
			this.Y(y);
			this.Z(z);
			this.W(w);
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
		
		public void W(double w) {
			this.mW = w;
			return;
		}

		public double W() {
			return this.mW;
		}
		

	}
	
	
	private static class Weight { // weight class

		private double mX = 0;
		private double mY = 0;
		private double mZ = 0;

		private Ssin ssin = null;
		
		public Weight(double x, double y) {
			this.X(x);
			this.Y(y);
		}

		public Weight(double x, double y, double z) {
			this.X(x);
			this.Y(y);
			this.Z(z);
		}
		public Weight(double x, double y, Ssin ssin) {
			this.X(x);
			this.Y(y);
			this.ssin = ssin;
		}
		
		public Weight(double x, double y, double z, Ssin ssin) {
			this.X(x);
			this.Y(y);
			this.Z(z);
			this.ssin = ssin;
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
		

	}
	
	private static class Data { // data class

		private double mX = 0;
		private double mY = 0;
		private int mCluster = 0;

		public Data() {
			
			return;
		}
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
	public static double max(double[] vec) {
		/**
		 * Function that finds the max val in
		 * a vector
		 */
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
