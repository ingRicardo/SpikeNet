/**
 * 
 */
package nnets;

import java.util.ArrayList;



/**
 * @author TIJUANA
 *
 */
public class ClusteringV1y {

	/**
	 * Using equation 5.9
	 * output clusters : 5
	 * Receptive Fields :1
	 * Type of Synapses: simples
	 * Training: weight and time adaptation
	 */
	private static ArrayList<Data> dataSet = new ArrayList<Data>();
	
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
	
	private static final int NUM_CLUSTERS = 5;
	private static final int TOTAL_DATA = SAMPLES.length;

	private static double [][] WEIGHTS = null;
	private static double [][] TAU = null;
	
	private static  double[][] nSAMPLES = new double [SAMPLES.length][SAMPLES[0].length];
	
	// Initial parameters
	
	private static int out_neu; //output neurons
	private static double rho; // time
	private static int conj; //training data
	private static int[] rf; // number of neurons by input
	private static double [] sigma; 
	private static double f_cut;
	private static int maxd;
	
	//Learning parameters
	private static double tau_max; // max time
	private static double tau_min; // min time
	private static double w_max; //max weight
	private static double w_min; // min weight
	private static int max_epoch; // codification interval
	private static int t_max; // max training time
	private static double eta_w; // learning rate by weights
	private static double beta ; // learning window of weights
	private static int nu; // neighborhood
	private static double dx;// weight motion
	
	// parameters 
	private static double kapa; // learning window
	private static int in_neu; // number of neurons by input
	private static double teta; // threshold
	
	private static double w; // weights and time initialization
	private static double tau; //
	
	private static int neu; 
	private static double inst;
	
	public ClusteringV1y(String title) {
	
		System.out.println(title);
		ini();
		training();
		test();
	}

	public static void ini() {

		// TODO Auto-generated constructor stub
	
		
		out_neu =5;
		rho = 0.1;
		conj = 5;
		rf =new int[] { 4,4};
		sigma = new double [] { 1/(1.5* (rf[0]-2)) ,  1/(1.5* (rf[0]-2))};
		f_cut = 0.9;
		maxd = 10;
		
		tau_max=5.8;
		tau_min = 0.1;
		w_max= 6.98;
		w_min =0;
		max_epoch =5;
		t_max = 10;
		eta_w = 0.035;
		beta=0.26;
		nu = 3;
		dx = 2.8;
		
		kapa = 1- Math.pow(nu, 2)/(2* Math.log(beta/(beta+1)));
		//in_neu = rf[0] + rf[1];
		in_neu = rf.length;
		teta=12;
		
		w = w_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) * (w_max-w_min);
		tau = tau_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) *(tau_max-tau_min); 
		
		WEIGHTS = new double[in_neu][TOTAL_DATA];
		TAU = new double[in_neu][TOTAL_DATA];
		
		// in_neu *  dataSet.size() = 8 * 10
		
		for (int i=0; i< in_neu; i++) {
			for (int j = 0; j < TOTAL_DATA; j++) {
				WEIGHTS[i][j] = (w_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) * (w_max-w_min));
				TAU[i][j] = (tau_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) *(tau_max-tau_min));
			}
			
		}
		
		
		System.out.println(" < weights > ");
	/*	for (int i=0; i< WEIGHTS.length; i++) {
			for (int j=0; j< WEIGHTS[0].length; j++) {
				System.out.println(" WEIGHTS[ "+i+" ][ "+j+" ] -> "+ WEIGHTS[i][j]);
				System.out.println( "TAU[ "+i+" ][ "+j+" ]  ->"+ TAU[i][j] );
			}
		}
		*/
		
		Data newData = null;
		int sampleNumber =0;
		double wx =0, wy=0;
		while( dataSet.size() < TOTAL_DATA) {
			wx = (w_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) * (w_max-w_min));
			wy = (w_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) * (w_max-w_min));
			newData = new Data(SAMPLES[sampleNumber][0], SAMPLES[sampleNumber][1]);
	
			dataSet.add(newData);
			sampleNumber++;
		}
		for (int i=0; i < dataSet.size(); i++) {
			System.out.println(dataSet.get(i).X() +" ,"+ dataSet.get(i).Y() );

		}
		System.out.println(" ");
		
		
		
	}
	public static void training() {
		System.out.println("training");
	
		kodieren_rf();// encoding data
		
		
		double [] tr = new double[nSAMPLES[0].length];
		int ctrl =0;
		double delta_t;
		while(ctrl <= max_epoch) {
			for(int i=0; i < conj; i++) {
				tr = new double[nSAMPLES[0].length];
				
				for (int j =0; j < nSAMPLES[0].length; j++) {
	
					tr[j] = nSAMPLES[i][j];
						
				}
				
				wer_spuckt(tr);
				
				System.out.println("neu -> "+ neu);
			//	dataSet.get(i).cluster(neu);
				System.out.println("isnt -> "+ inst);
		//		System.out.println("tr length -> " + tr.length );
				for (int x =0 ; x < in_neu; x++) {
					
					delta_t= tr[x] - inst;
				//	System.out.println("delta_t ->  " + delta_t);
					
					WEIGHTS[x][neu] = WEIGHTS[x][neu] + eta_w * ((1-beta)* 
					Math.exp( -1* ( Math.pow((delta_t+dx), 2) ) / (2*kapa-2) )-beta);
					
					TAU[x][neu] = tau_max * (0.2452* Math.exp(-0.07922 *(WEIGHTS[x][neu]/w_max)  ) 
							+ 0.09115 * Math.exp(  2.149*(WEIGHTS[x][neu]/w_max)   )  );
					
					System.out.println("WEIGHTS[ "+x+" ][ "+neu+" ] -> "+ WEIGHTS[x][neu]);
					System.out.println("TAU[ "+x+" ][ "+neu+" ] -> "+ TAU[x][neu]);
					
					if ( WEIGHTS[x][neu] > w_max)
						 WEIGHTS[x][neu] = w_max;
					
					if ( WEIGHTS[x][neu] < w_min)
						WEIGHTS[x][neu] = w_min;
					
					if ( TAU[x][neu] > tau_max)
						 TAU[x][neu] = tau_max;
					
					if ( TAU[x][neu] < tau_min)
						TAU[x][neu] = tau_min;
				}
				
			}
			System.out.println("ctrl -> " + ctrl);
			ctrl++;
			
			
			teta = teta + ((0.3*teta)/max_epoch);
		}//while
		
		
	}
	
	public static void test() {
		// testing
				System.out.println(" <-testing -> ");
				double [] tr = new double[nSAMPLES[0].length];
			
				for(int i=0; i <NUM_CLUSTERS ; i++) {
					tr = new double[nSAMPLES[0].length];
					
					for (int j =0; j < nSAMPLES[0].length; j++) {

						tr[j] = nSAMPLES[i][j];
							
					}
					
					wer_spuckt(tr);
					
					System.out.println("neu -> "+ neu);
					dataSet.get(i).cluster(neu);
					
				}
				for(int i=0; i< dataSet.size(); i++) {
					System.out.println("clus -> "+ dataSet.get(i).cluster());
					
				}
				
	}
	public static void wer_spuckt(double [] tr) {
		System.out.println(" <- wer_spuckt -> ");
		
		double t =0;
		neu = 0;
		double [] out = null;
		double  sai=0;
		double dt =0;

		while(neu == 0 &&  t <= t_max) {
			out = new double[ TOTAL_DATA];		
			for(int j=0; j <TOTAL_DATA; j++) {
				//	System.out.println("");
				for (int i=0 ; i  < tr.length; i++) {
					
						dt = ((t - tr[i]) / TAU[i][j]);
						sai = WEIGHTS[i][j]* dt* Math.exp(1- dt);
											
						if (sai < 0)
							sai=0;					
						
						out[j] = out[j] + sai;
				}
				
			}
			
			double max= max(out);
			if (max >= teta ) {
				
				for(int i =0; i< out.length; i++) {
					
					if (max == out[i]) {
						neu =i;
						inst = t;
					
					}
				}
				
			}
				
			t=t+rho;
		}//while
	}
	
	public static void kodieren_rf() {
		System.out.println("<-Encoding Data->");
		double [] ct = {rf[0], rf[1]};		
		double [] tmp = new double[SAMPLES.length]; 
		double [] tmpX = new double[tmp.length];
		double [] tmpY = new double[tmp.length];
		
		if (ct[0] > 3 || ct[1] > 3 ) {
			for (int j =0; j < SAMPLES[0].length; j++) {
				tmp = new double[SAMPLES.length]; 
				double range =0;
				for (int i=0; i < SAMPLES.length; i++) {
					tmp[i] = SAMPLES[i][j];
				}

				range = max(tmp) - min(tmp);
				for (int i=0; i < SAMPLES.length; i++) {
					
					if (j ==0) {
						tmpX[i] = (tmp[i] - min(tmp))/range;
					}
					if (j==1) {
						tmpY[i] = (tmp[i] - min(tmp))/range;
					}
				}
		  }
			
		Data newData = null;
		int sampleNumber =0;
		dataSet = new ArrayList<Data>();
		while( dataSet.size() < TOTAL_DATA) {
			newData = new Data(tmpX[sampleNumber], tmpY[sampleNumber]);
			dataSet.add(newData);
			sampleNumber++;
		}
		for (int i=0; i < dataSet.size(); i++) {
			System.out.println(dataSet.get(i).X() +" ,"+ dataSet.get(i).Y());
			nSAMPLES[i][0] = dataSet.get(i).X();
			nSAMPLES[i][1] = dataSet.get(i).Y();
		}
	}//if
		
		
	}
	
	
	private static class Data { // data class

		private double mX = 0;
		private double mY = 0;
		private double wX =0;
		private double wY =0;
		
		private int mCluster = 0;

		public Data() {
			
			return;
		}
		public Data(double x, double y) {
			this.X(x);
			this.Y(y);
		}

		public Data(double x, double y, double wx, double wy) {
			this.X(x);
			this.Y(y);
			this.wX(wx);
			this.wY(wy);
		}
		
		public void wX(double wx) {
			this.wX = wx;
			return;
		}

		public double wX() {
			return this.wX;
		}

		public void wY(double wy) {
			this.wY = wy;
			return;
		}

		public double wY() {
			return this.wY;
		}

		
		public void Y(double y) {
			this.mY = y;
			return;
		}

		public double Y() {
			return this.mY;
		}

		
		public void X(double x) {
			this.mX = x;
			return;
		}

		public double X() {
			return this.mX;
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
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ClusteringV1y cluster = new ClusteringV1y("Cluster SNN");
	}

}
