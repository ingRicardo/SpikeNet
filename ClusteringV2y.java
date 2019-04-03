/**
 * 
 */
package nnets;

import java.util.ArrayList;



/**
 * @author TIJUANA
 *
 */
public class ClusteringV2y {

	/**
	 * Using equation 5.9
	 * output clusters : 5
	 * Receptive Fields :1
	 * Type of Synapses: simples
	 * Training: weight and time adaptation
	 */
	private static ArrayList<Data> dataSet = new ArrayList<Data>();
	
	private static final double[][] SAMPLES = { 
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 },
			{ (Math.random()*((20.0-1.0)+1))+ 1.0, (Math.random()*((20.0-1.0)+1))+ 1.0 } };
	
	private static final int NUM_CLUSTERS = 5;
	private static final int TOTAL_DATA = SAMPLES.length;

	private static double [] WEIGHTS = null;
	private static double [] TAU = null;
	
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
	
	public ClusteringV2y(String title) {
	
		System.out.println(title);
		ini();
		training();
		test();
		print();
	}

	public static double doubRandomW() {
		
		 // define the range 
        int max = in_neu; 
        int min = out_neu; 
        double range = max - min + 1; 
        //(w_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) * (w_max-w_min));
            double rand = w_min + (((Math.random() * range) + min)+ out_neu)* (w_max-w_min); 
  
            return rand;
     
	}
	
	public static double doubRandomT() {
		
		 // define the range 
       int max = in_neu; 
       int min = out_neu; 
       double range = max - min + 1; 
       //(tau_min + ( (Math.random()*((in_neu-out_neu)+1))+ out_neu) *(tau_max-tau_min));
           double rand = tau_min + (((Math.random() * range) + min)+ out_neu)*(tau_max-tau_min); 
 
           return rand;
    
	}
	
	public static void ini() {

		// TODO Auto-generated constructor stub
	
		
		out_neu =5;
		rho = 0.2;
		conj = 5;
		rf =new int[] { 5,5};
		sigma = new double [] { 1/(1.5* (rf[0]-2)) ,  1/(1.5* (rf[0]-2))};
		f_cut = 0.9;
		maxd = 10;
		
		tau_max=5.8;
		tau_min = 0.1;
		w_max= 6.98;
		w_min =0;
		max_epoch =10;
		t_max = 10;
		eta_w = 0.665;
		beta=0.56;
		nu = 9;
		dx = 4.8;
		
		kapa = 1- Math.pow(nu, 2)/(2* Math.log(beta/(beta+1)));
		in_neu = rf[0] + rf[1];
		//in_neu = rf.length;
		teta=4;
		
		WEIGHTS = new double[in_neu];
		TAU = new double[in_neu];
				
		for (int i=0; i< in_neu; i++) {
			
				WEIGHTS[i]= doubRandomW();
				TAU[i] = doubRandomT();
		}
		
		
		System.out.println(" < weights > ");
		for (int i=0; i< WEIGHTS.length; i++) {			
				System.out.println(" WEIGHTS[ "+i+ "] -> "+ WEIGHTS[i]);
				System.out.println( "TAU[ "+i+" ]  ->"+ TAU[i] );			
		}
		
		
		Data newData = null;
		int sampleNumber =0;
		while( dataSet.size() < TOTAL_DATA) {

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
			System.out.println("ctrl -> " + ctrl);
			for(int i=0; i < conj; i++) {
				tr = new double[nSAMPLES[0].length];
				
				for (int j =0; j < nSAMPLES[0].length; j++) {
	
					tr[j] = nSAMPLES[i][j];
						
				}
			
				wer_spuckt(tr);
				
				System.out.println("neu -> "+ neu);

				System.out.println("isnt -> "+ inst);
					
					delta_t= (tr[0]+tr[1]) - inst;
					
					WEIGHTS[neu] = WEIGHTS[neu] + eta_w * ((1-beta)* 
					Math.exp( -1* ( Math.pow((delta_t+dx), 2) ) / (2*kapa-2) )-beta);
					
					TAU[neu] = tau_max * (0.2452* Math.exp(-0.07922 *(WEIGHTS[neu]/w_max)  ) 
							+ 0.09115 * Math.exp(  2.149*(WEIGHTS[neu]/w_max)   )  );
					
					System.out.println("WEIGHTS[ "+neu+" ] -> "+ WEIGHTS[neu]);
					System.out.println("TAU[ "+neu+" ] -> "+ TAU[neu]);
					
					if ( WEIGHTS[neu] > w_max)
						 WEIGHTS[neu] = w_max;
					
					if ( WEIGHTS[neu] < w_min)
						WEIGHTS[neu] = w_min;
					
					if ( TAU[neu] > tau_max)
						 TAU[neu] = tau_max;
					
					if ( TAU[neu] < tau_min)
						TAU[neu] = tau_min; 
				
			}
			
			ctrl++;
			
			
			teta = teta + ((0.3*teta)/max_epoch);
		}//while
		
		
	}
	
	public static void test() {
		// testing
				System.out.println(" <-testing -> ");
				
				System.out.println(" < weights > ");
				for (int i=0; i< WEIGHTS.length; i++) {			
						System.out.println(" WEIGHTS[ "+i+ "] -> "+ WEIGHTS[i]);
						System.out.println( "TAU[ "+i+" ]  ->"+ TAU[i] );			
				}
				
				
				
				
				double [] tr = new double[nSAMPLES[0].length];
		
				for(int i=0; i <TOTAL_DATA ; i++) {
					tr = new double[nSAMPLES[0].length];
					
					for (int j =0; j < nSAMPLES[0].length; j++) {

						tr[j] = nSAMPLES[i][j];

					}
					
					wer_spuckt(tr);
					
					System.out.println("neu -> "+ neu);
					if(neu > NUM_CLUSTERS)
						neu = 0;
					
					if(neu == NUM_CLUSTERS)
						neu = NUM_CLUSTERS-1;
					
					dataSet.get(i).cluster(neu);
					
					
				}
				for(int i=0; i< dataSet.size(); i++) {
					System.out.println("clus -> "+ dataSet.get(i).cluster());
					
				}
				
	}
	
	public static void print() {
		
		for (int i = 0; i < NUM_CLUSTERS; i++) {

			System.out.println("Cluster " + i + " includes: ");

			for (int j = 0; j < TOTAL_DATA; j++) {

				if (dataSet.get(j).cluster() == i) {
					System.out.println(" (" + dataSet.get(j).X() + ", " + dataSet.get(j).Y() + ")");
				}
			}
			System.out.println();

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
						for(int p=0; p <in_neu; p++) {
							double sum=0;
							for (int i=0 ; i  < tr.length; i++) {
								dt = ((t - tr[i]) / TAU[p]);
								sai = WEIGHTS[p]* dt* Math.exp(1- dt);
													
								if (sai < 0)
									sai=0;					
								
								sum = sum + sai;
								
							}
							out[p] = sum;

					}//data
					
					double max= max(out);
				
					if (max >= teta ) {
						System.out.println("max -> "+ max + " teta -> " + teta);
						for(int j =0; j< out.length; j++) {
							System.out.println("out[ "+j+" ] - > " + out[j]);
							if (max == out[j]) {
								neu =j;
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
		ClusteringV2y cluster = new ClusteringV2y("Cluster SNN");
	}

}
