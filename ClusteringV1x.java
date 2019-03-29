/**
 * 
 */
package nnets;

import java.util.ArrayList;

/**
 * @author TIJUANA
 *
 */
public class ClusteringV1x {

	/**
	 * @param args
	 */
	
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
	private static int out_neu ; // output neurons
	private static double rho; //time
	private static int conj; //training data
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
	private static double dx;//motion of learning window
	
	//----------------
	private static double kapa; // number for learning window
	private static double in_neu; // number of neurons by input
	private static int ssin ; // number of sub-synapses
	private static double teta;// threshold
	
	//----------------
	private static int typ; // curve center
	private static double inst;
	private static int neu;
//	private static double []
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		initialize();
		Training();
	}
	
	public static void initialize() {
		out_neu = 5;
		rho = 0.1;
		conj = SAMPLES.length/2; //training data
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
		dx =2.3;
		//-------------------------
		kapa = 1 - Math.pow(nu, 2) / (2 * Math.log(beta / (beta+1)));
		in_neu = rf[0] + rf[1];
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
		
		for (int i=0; i <dataSet.size(); i++) {
			
			System.out.println(dataSet.get(i).X() + " , "+ dataSet.get(i).Y());
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
		
		System.out.println(" ");
		for (int i=0; i <weightSet.size(); i++) {
			
			System.out.println(weightSet.get(i).X() + " , "+ weightSet.get(i).Y());
			System.out.println(" ssin -> "+  weightSet.get(i).ssin.X());
			System.out.println(" ssin -> "+  weightSet.get(i).ssin.Y());
			System.out.println(" ssin -> "+  weightSet.get(i).ssin.Z());
			System.out.println(" ssin -> "+  weightSet.get(i).ssin.W());
		}
	}
	
	public static void Training() {
		
		kodieren_rf(); // encode input
		
		//------
	/*	for (int z =0; z< DelayDataSet.size(); z++) {
			
			[neu,inst]=wer_schiesst(tr,t_max,w,d,tau,rho,teta,0); % quem disparou?
		}*/
	int ctrl_1 =0;
	while ( ctrl_1 <= max_epoch) {
		wer_schiesst();
		double [] delta_t = new double[ssin]; 
		if (neu !=0) {
			for (int i =0; i< DelayDataSet.size(); i++) {
				
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
						for (int p =0; p < ssin; p++) {
							delta_t[p] = (DelayDataSet.get(i).X() +  DelayDataSet.get(i).Y()) + d[p] -inst;
						}
					/*	weightSet.get(neu).X(
											weightSet.get(neu).X()+ eta * 
												( (1 +beta) *
												  Math.exp(-1 * (Math.pow((delta_t[0]+2.3), 2) / (2*kapa-2) ))-beta));
						
						weightSet.get(neu).Y(
								weightSet.get(neu).Y()+ eta * 
									( (1 +beta) *
									  Math.exp(-1 * (Math.pow((delta_t[0]+2.3), 2) / (2*kapa-2) ))-beta));
						*/
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
						
						
						
					}
					
					
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
					
				}
				
				
				
			}
			
		}
		System.out.println("ctrl_1 -> "+ ctrl_1);
		ctrl_1++;
	 }//while
	}
	
	public static void wer_schiesst() {
		int [] indX = new int[DelayDataSet.size()];
		int [] indY = new int[DelayDataSet.size()];
		for (int z =0; z< DelayDataSet.size(); z++) {
			
			if ( DelayDataSet.get(z).X() != -1 ) {
				System.out.println(z);
				indX[z] = z;
			}else {
				indX[z] = 100;
			}

			if ( DelayDataSet.get(z).Y() != -1 ) {
				System.out.println(z);
				indY[z] = z;
			}else {
				indY[z] = 100;
			}

			
		}
		
		System.out.println( "\n X index fired ");
		for (int i =0; i < indX.length; i++)
			System.out.println(indX[i]);
		
		System.out.println( "\n Y index fired ");
		for (int i =0; i < indY.length; i++)
			System.out.println(indY[i]);

			double t =0; 
		
	//	for ( int j =0; j< SAMPLES[0].length; j++) { // aus
			double dt=0, dtY=0;
			double dtt = 0, dttY = 0, dttZ = 0, dttW = 0;
			double sai=0;
			
			neu =0;
			while (neu == 0 && t <= t_max ) {
			
			double [] out = new double[DelayDataSet.size()];
			int i =0;
			for ( i =0; i < DelayDataSet.size(); i++) { //tam
												
				if (indX[i] == i) {
					
					dt = t - (DelayDataSet.get(i).X() + DelayDataSet.get(i).Y());
					
					for (int k = 0; k < ssin; k++ ) {
																
						dtt = (dt - d[k]) / tau;
																		
						sai = 	weightSet.get(i).X() * dtt * Math.exp(1 - dtt) +
								weightSet.get(i).Y() * dtt * Math.exp(1 - dtt) +						
								weightSet.get(i).ssin.X() * dtt * Math.exp(1 - dtt) +
								weightSet.get(i).ssin.Y() * dtt * Math.exp(1 - dtt) +
								weightSet.get(i).ssin.Z() * dtt * Math.exp(1 - dtt) +
								weightSet.get(i).ssin.W()* dtt * Math.exp(1 - dtt);
				
						System.out.println("sai -> "+ sai);
						if (sai < 0)
							sai =0;
						
						out[i] = out[i] +sai;	
						System.out.println("out[ "+i+" ] -> "+ out[i] );
					}//ssin
						
				}//indx
			}// tam
			System.out.println( " out lenght -> "+ out.length +" max "+ max(out) + " teta -> "+ teta);
			if (max(out) >= teta) {
				System.out.println("i ---> "+i);
				for (int j=0; j< out.length; j++) {
					
					if (max(out) == out[j]) {
						neu = j;
						inst = t;
						System.out.println("  neu ------> "+ neu + " inst --> " + inst);
					}
						
				}
				
			}
			t = t + rho;
			}
			
			
	//	}
		
	}
	
	
	public static void kodieren_rf() {
		System.out.println("kodieren_rf() \n");
		int t_neur = rf[0] + rf[1];
		if ( rf[0] >= 3 && rf[1] >= 3) {
			
			double [][] aus = new double[SAMPLES.length][SAMPLES[0].length];
			int [] ct = {rf[0], rf[1]};
			
			double [] tmpVec = new double[SAMPLES[0].length];
			double [] tmpVecX = new double[SAMPLES[0].length];
			double [] tmpVecY = new double[SAMPLES[0].length];
			double range, min, max;
		
			int j =0;
			for ( j =0; j< SAMPLES[0].length; j++) {
				tmpVec = new double[SAMPLES.length];
				System.out.println(" ");
				for (int i =0; i < SAMPLES.length; i++) {
					tmpVec[i] =SAMPLES[i][j];
					System.out.println(" tmpVec[ "+i+" ] -> "+ tmpVec[i]);
				}
				min = min(tmpVec);
				max = max(tmpVec);
				range= max - min;
				System.out.println("range -> "+ range);
				
				for (int i =0; i < tmpVec.length; i++) {
					tmpVec[i]=(tmpVec[i] - min)/range;
				}
				
				System.out.println(" <- new vals -> \n");
				
				if ( j == 0)
					tmpVecX = new double[SAMPLES.length];
				
				if (j == 1)
					tmpVecY = new double[SAMPLES.length];
				
				for (int i =0; i < tmpVec.length; i++) {
					
					if ( j == 0)
						tmpVecX[i] = 	tmpVec[i];
					
					if ( j == 1)
						tmpVecY[i] = 	tmpVec[i];
					
					
					System.out.println(" tmpVec[ "+i+" ] -> "+ tmpVec[i]);
										
					
				}
				
			}
			Data newData = null;
			int numberSample =0;
			while(NormalDataSet.size() < TOTAL_DATA) {
				
				newData = new Data(tmpVecX[numberSample], tmpVecY[numberSample]);
				NormalDataSet.add(newData);
				numberSample++;
			}
			
			System.out.println("\n< - NormalDataSet ->");
			
			for (int i =0; i < NormalDataSet.size(); i++) {
				System.out.println("");
				System.out.println(" X " + NormalDataSet.get(i).X() + " Y " + NormalDataSet.get(i).Y());
			}
		double cdis =0;  // distance between centers
		if (typ == 1) {
			
			cdis = 1/rf[0];
			
		} else if ( typ == 0) {
			
			cdis = 1/(rf[0]-2);
		}
		
		for (int i =0; i< rf.length; i++) {
			
			ct[i] = (int) ((i -0.5) * cdis);
		}
		
		// getting delays 
		double auxX =0, auxY =0;

		double [] ausX = new double[TOTAL_DATA] ;
		double [] ausY = new double[TOTAL_DATA] ;

		System.out.println("\n<-- getting delays --> ");
			Data newDelay = null;
			
		for (int i =0; i < NormalDataSet.size(); i++) {
			System.out.println("");

					 auxX = maxd - maxd *
							 Math.exp(-1 * (Math.pow((NormalDataSet.get(i).X() -ct[0]), 2)
									 / (2* Math.pow(sigma[0], 2))));
 
					 auxY = maxd - maxd *
							 Math.exp(-1 * (Math.pow((NormalDataSet.get(i).Y() -ct[1]), 2)
									 / (2* Math.pow(sigma[1], 2))));

				
					 ausX[i] = Math.round(auxX*(1/rho)) * rho; // round
					 ausY[i] = Math.round(auxY*(1/rho)) * rho; // round
				 

				 System.out.println(" ausX[ "+i+" ] ->  " +ausX[i]  +" ausY[ "+i+" ] ->  " +ausY[i]  + " maxd * f_cut -> " + (maxd * f_cut));
				 
				 if ((ausX[i] > ( maxd * f_cut)) && (ausY[i] > ( maxd * f_cut))) { // firing or not
					 System.out.println("** firing **");
					 ausX[i] = -1;			// delay of -1
					 ausY[i] = -1;
				 }
				 
			
			newDelay = new Data(ausX[i], ausY[i]);		
			DelayDataSet.add(newDelay);
	
			
		}
		 
		System.out.println(" <- ausX values -> ");
		
		for (int i =0; i< ausX.length; i++)
			System.out.println(ausX [i]);
		
		System.out.println(" <- ausY values -> ");
		
		for (int i =0; i< ausY.length; i++)
			System.out.println(ausY [i]);

		
		System.out.println("\n< - DelayDataSet ->");
		
		for (int i =0; i < DelayDataSet.size(); i++) {
			System.out.println("");
			System.out.println(" X " + DelayDataSet.get(i).X()+ " Y " + DelayDataSet.get(i).Y());
		}
		
	
		}else {
			System.out.println("Incorrect size of receptive fields it must be greater than 2");
		}
		

	}
	
	
	
	private static class Ssin {

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
	
	
	private static class Weight {

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
	
	private static class Data {

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
