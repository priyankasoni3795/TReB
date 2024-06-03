
package METO_Project;

import java.lang.*;
import java.util.Scanner;

class Criteria_fog{
	public int quota;
    public double latency;
    public double energyConsumption;
    public Criteria_fog(double x,double y){
        latency=x;
        energyConsumption=y;
    }
}

class Criteria_IoT{
    public double energy;
    public double deadline;
    public Criteria_IoT(double x,double y){
    	energy=x;
        deadline=y;
    }
}

class PAIR{
    public double val;
    public int no;
    public PAIR(double x,int y){
        val=x;
        no=y;
    }
}


public class METO {
	static Criteria_IoT[] fogDevices;
	static Criteria_fog[] tasks;
    static double[][][] decisionMatrixTask;
    static double[][][] decisionMatrixFN;
    static double[][] weightMatrixTask;
    static double[][] weightMatrixFN;
    static double[][] performanceTask;
    static double[][] performanceFN;
    static Integer[] Q;
    static Integer[] numberOfDependentDevices;
    static int noOfFN;
	static int noOfTask;
	
	public static void main(String[] args) {
        int n=250;		//No of devices
        int m=5;		//No of fog nodes
        Scanner sc = new Scanner (System.in);
        //taking all required inputs
        System.out.println ("please provide the following inputs::\n");
//        System.out.println ("enter the computational demand of the task in million cycle::");
//        double computational_demand = sc.nextDouble();
        System.out.println ("enter frequency in megahertz::");
        double frequency = 16; //sc.nextDouble ();
        System.out.println ("enter channel frequency(bandwidth) in megahetz::");
        double bandwidth = 10; //sc.nextDouble ();
        System.out.println ("enter Distance of IOT device to fog node in meter::");
        int distance_to_fog_node = 100; //sc.nextInt ();
        System.out.println ("enter the power of iot device in watts::");
        double power_of_iot_device = 0.5; //sc.nextDouble ();
        System.out.println ("enter the power of fog node in watts :: ");
        double power_of_fogDevice = 1.5; //sc.nextDouble ();
        System.out.println ("enter the output size in Kb::");
        double output_size = 20; //sc.nextDouble ();
        System.out.println ("enter the computational power of fog node::");
        double computational_power_of_fognode = 0.35; //sc.nextDouble ();
        
        
        
        
        noOfTask=n;
        noOfFN=m;
        tasks=new Criteria_fog[n];
	    fogDevices=new Criteria_IoT[m];
	    double x;
        double y;
        
        Q=new Integer[noOfFN];
        
        for(int i=0;i<m;i++){
            x=calculate_energy_fog(distance_to_fog_node, power_of_fogDevice, power_of_iot_device, bandwidth, output_size) ;
            System.out.println ("deadline of fognode " + (i+1) + " : ");
    		y = sc.nextDouble ();
    		
    		// Q -- Quota of each FN
    		System.out.println ("Quota of fognode " + (i+1) + " : ");
    	    Q[i]=sc.nextInt(); // store some value
            
    	    Criteria_IoT item = new Criteria_IoT(x ,y);
            fogDevices[i]=item;
        }
        
        numberOfDependentDevices = new Integer[noOfTask];
        for(int i=0; i<n; i++) {
        	numberOfDependentDevices[i] = 0;
        }
        
        double input_for_tasks[][] = new double[noOfTask][2];
        for(int i=0;i<4;i++){
        	System.out.println ("enter computational demand of  task " + (i+1) +" :: " );
        	input_for_tasks[i][0] = sc.nextDouble ();
        	System.out.println ("enter inputSize of task " + (i+1) + " :: " );
        	input_for_tasks[i][1] = sc.nextDouble ();
        	System.out.println ("enter number Of Dependent Devices for task " + (i+1) + " :: " );
        	numberOfDependentDevices[i] = sc.nextInt ();
        }
        for(int i=4;i<n;i++){
        	System.out.println ("enter computational demand of  task " + (i+1) +" :: " );
        	input_for_tasks[i][0] = sc.nextDouble ();
        	System.out.println ("enter inputSize of task " + (i+1) + " :: " );
        	input_for_tasks[i][1] = sc.nextDouble ();
        }
        
        for(int i=0;i<n;i++) {
        	double temp[] = calculate_delay_and_energy(distance_to_fog_node, power_of_iot_device, bandwidth, output_size, power_of_fogDevice, i, input_for_tasks);
            x= temp[0];
            y= temp[1];
            Criteria_fog item = new Criteria_fog(x ,y);
            tasks[i]=item;
        }
        
		decisionMatrixTask=new double[n][m][2];
		decisionMatrixFN=new double[m][n][2];
		
		weightMatrixTask=new double[n][2];
		weightMatrixFN=new double[m][2];
		
		performanceTask=new double[n][m];
		performanceFN=new double[m][n];
		
		createDecisionMatrix();
//		print1();
		
		CRITIC();
//		print2();
		
		TOPSIS();
//		print3();
		
		matching();
		
		
		//total energy of the system
		double total_energy = 0.0;
		for(int i=0;i<n;i++){
	    	total_energy += tasks[i].energyConsumption;
	    }
	    for(int i=0; i<m; i++) {
	    	total_energy += fogDevices[i].energy;
	    }
	    System.out.println("\nTotal Energy for " + noOfTask +" devices = "+ total_energy);
	    
	    double total_latency = 0.0;
	    for(int i=0; i<n; i++) {
	    	total_latency += tasks[i].latency;
	    }
	    System.out.println("\nTotal Latency for " + noOfTask +" devices = "+ total_latency);

	}
	
	
	public static double log10 (int x)
	  {
	    return Math.log (x) / Math.log (10);
	  }
	
	public static double log2 (double x)
	  {
	    return Math.log (x) / Math.log (2);
	  }
	
	public static double[] calculate_delay_and_energy(int distance_to_fog_node, double power_of_iot_device, double bandwidth, double output_size, double power_of_fogDevice, int i, double t[][]) {
		double computational_demand  = t[i][0];
		double input_size = t[i][1];
        
        double log10distance_to_fog_node = log10 (distance_to_fog_node);

        double PLIJ = 38.02 + 20 * log10distance_to_fog_node;
        //   System.out.println(PLIJ);
        double power_to_calculate = -(PLIJ / 10);
        double channel_gain = Math.pow (10, power_to_calculate);

        double noise_power = Math.pow (10, -10);

        double in_the_log = 1 + (power_of_iot_device * channel_gain) / noise_power;

        double log2in_the_log = log2 (in_the_log);
        double uplink_data_rate = bandwidth * log2in_the_log;								//U(i,j)

        double Transmission_delay = input_size / uplink_data_rate;
        
        double receiving_delay = output_size / uplink_data_rate;
        
        double fog_node_capacity = 6 * Math.pow(10, 9);
        double temp_quota = 50.0;
        double capacity_of_VRU = fog_node_capacity/temp_quota ;
        
        double computational_delay = (computational_demand) / (capacity_of_VRU);
        
        double total_delay = Transmission_delay + receiving_delay + computational_delay;
        
        double energy_IoT = power_of_iot_device*(Transmission_delay + receiving_delay);
        		
        
      
        double res[] = new double[2];
        res[0] = total_delay;
        res[1] = energy_IoT;
        return res;
	}
	
	public static double calculate_energy_fog(int distance_to_fog_node, double power_of_fogDevice, double power_of_iot_device, double bandwidth, double output_size) {
		double computational_demand = 300;
    	double input_size = 500;
        
        double log10distance_to_fog_node = log10 (distance_to_fog_node);

        double PLIJ = 38.02 + 20 * log10distance_to_fog_node;
        //   System.out.println(PLIJ);
        double power_to_calculate = -(PLIJ / 10);
        double channel_gain = Math.pow (10, power_to_calculate);

        double noise_power = Math.pow (10, -10);

        double in_the_log = 1 + (power_of_iot_device * channel_gain) / noise_power;

        double log2in_the_log = log2 (in_the_log);
        double uplink_data_rate = bandwidth * log2in_the_log;								//U(i,j)

        double Transmission_delay = input_size / uplink_data_rate;
        
        double receiving_delay = output_size / uplink_data_rate;
        
        double fog_node_capacity = 6 * Math.pow(10, 9);
        double temp_quota = 50.0;
        double capacity_of_VRU = fog_node_capacity/temp_quota ;
        
        double computational_delay = (computational_demand) / (capacity_of_VRU);

        double energy_fog =	power_of_fogDevice * (Transmission_delay + receiving_delay) + power_of_iot_device * (computational_delay);
		return energy_fog;
	}
	
	private static void createDecisionMatrix(){
		

		// creating decision matrix for Tasks of IOT devices
		// matrix dTask size = noOfFN x 2
		
		for(int k=0;k<noOfTask;k++){  
		    double[][] dTask=new double[noOfFN][2];  
			for(int i=0;i<noOfFN;i++){
				
				dTask[i][0]=fogDevices[i].energy;		//+Math.abs(Math.random());
				dTask[i][1]=fogDevices[i].deadline;		//+Math.abs(Math.random());
			}
			decisionMatrixTask[k]=dTask;
		}
		

		// creating decision matrix for FN
		// matrix dFN size = noOfTask x 2
		for(int k=0;k<noOfFN;k++){
			
			double[][] dFN=new double[noOfTask][2];  
			for(int i=0;i<noOfTask;i++){
				
				dFN[i][0]=tasks[i].latency;					//+Math.abs(Math.random());
				dFN[i][1]=tasks[i].energyConsumption;		//+Math.abs(Math.random());
			}
			decisionMatrixFN[k]=dFN;
		}

	}
	
    private static double max(double x,double y){
	    if(x>y) return x;
	    else return y;
	}
	private static double min(double x,double y){
	    if(x<y) return x;
	    else return y;
	}



	private static void CRITIC(){

    	

    	for(int l=0;l<noOfTask;l++){
        	
        	int n=noOfFN;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixTask[l];
        	// Normalize the decision matrix of an agent ‘a’ as per Eq. (16).
        
        	double[] best=new double[c];
        	double[] worst=new double[c];
        	for(int j=0;j<c;j++){
            	best[j]=B[0][j];
            	worst[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	best[j]=max(best[j],B[i][j]);
                	worst[j]=min(worst[j],B[i][j]);
            	}
        	}
        	
        	for(int j=0;j<c;j++){
   		 		for(int i=0;i<n;i++){
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-worst[j]);  // Executing Eq 16
            	}
        	}
            
            

        	//Evaluate the standard deviation σk for each criterion in the normalized decision matrix
        	double[] SD=new double[c];
        	for(int j=0;j<c;j++){
            	
            	double[] v=new double[n];
            	for(int i=0;i<n;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v,n);
        	}
            

        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        
        	double[][] S=new double[c][c];
        	for(int i=0;i<c;i++){
        		
            	for(int j=0;j<c;j++){
            		double[] v1=new double[n];
            		double[] v2=new double[n];
                	for(int x=0;x<n;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<n;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2,n);
            	}
        	}
         

           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ‘a’
        	
        	double[] w=new double[c];
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<c;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
	     	
	     
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)
            	

            //Add Wa as next row in W.
   			weightMatrixTask[l]=w;     	        
    	}
    	
    	
// --------------------------------------------------------------------    	
    	
    	
    	for(int l=0;l<noOfFN;l++){
        	int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixFN[l];
        	// Normalize the decision matrix of an agent ‘a’ as per Eq. (16).
        	
        	double[] best=new double[c];
        	double[] worst=new double[c];
        	for(int j=0;j<c;j++){
            	best[j]=B[0][j];
            	worst[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	best[j]=max(best[j],B[i][j]);
                	worst[j]=min(worst[j],B[i][j]);
            	}
        	}
        	for(int j=0;j<c;j++){
   		 		for(int i=0;i<n;i++){
                	B[i][j]=(B[i][j]-worst[j])/(best[j]-worst[j]);  // Executing Eq 16
            	}
        	}


        	//Evaluate the standard deviation σk for each criterion in the normalized decision matrix
        	double[] SD=new double[c];
        	for(int j=0;j<c;j++){
            	
            	double[] v=new double[n];
            	for(int i=0;i<n;i++) v[i]=B[i][j];
            	SD[j]=calculateSD(v,n);
        	}


        	// constructing a criteria correlation matrix which is a symmetric matrix of size  c*c
        	
        
        	double[][] S=new double[c][c];
        	for(int i=0;i<c;i++){
        		
            	for(int j=0;j<c;j++){
            		double[] v1=new double[n];
            		double[] v2=new double[n];
                	for(int x=0;x<n;x++) v1[x]=(B[x][i]);
                	for(int x=0;x<n;x++) v2[x]=(B[x][j]);
                	S[i][j]=find_coefficient(v1,v2,n);
            	}
        	}


           // Determine each criterion weight wk calculated as per Eq. (18) and form the criteria weight vector Wa for agent ‘a’
        	
        	double[] w=new double[c];
        	float total_sum=0;
        	for(int i=0;i<c;i++){
            	float sum=0;
            	for(int j=0;j<c;j++)
                	sum+=(1-S[i][j]);
            	w[i]=SD[i]*sum;        // Executing Eq (17)
	     	}
			for(int i=0;i<c;i++)
				total_sum+=w[i];

        	for(int i=0;i<c;i++)
            	w[i]=(w[i]/total_sum);   //  Executing Eq (18)

            //Add Wa as next row in W.
   			weightMatrixFN[l]=w;     	        
    	}
   	}

	private static double find_coefficient(double[] X, double[] Y, int n){
   		double sum_X = 0, sum_Y = 0, sum_XY = 0;
   		double squareSum_X = 0, squareSum_Y = 0;
   		for (int i = 0; i < n; i++){
      		sum_X = sum_X + X[i];
      		sum_Y = sum_Y + Y[i];
      		sum_XY = sum_XY + X[i] * Y[i];
      		squareSum_X = squareSum_X + X[i] * X[i];
      		squareSum_Y = squareSum_Y + Y[i] * Y[i];
   		}
   		double corr = (float)(n * sum_XY - sum_X * sum_Y) / Math.sqrt((n * squareSum_X - sum_X * sum_X) * (n * squareSum_Y - sum_Y * sum_Y));
   		return corr;
	}

	private static double calculateSD(double[] data,int n) {
  		double sum = 0.0, mean, standardDeviation = 0.0;
  		int i;

  		for(i = 0; i < n; ++i) {
    		sum += data[i];
  		}

  		mean = sum / n;

  		for(i = 0; i < n; ++i) {
    		standardDeviation += Math.pow(data[i] - mean, 2);
  		}

  		return Math.sqrt(standardDeviation / n);
	}

	private static void TOPSIS(){
    	

    	for(int l=0;l<noOfTask;l++){
       		int n=noOfFN;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixTask[l];
            
    
                
                
        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}
        	
        	
        	

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixTask[l][j];  // Eq. (20)
        	}



           
                
                
        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}
            
           
                
        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}
        	
        	
           
                
           
                
                
        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}
            
//             System.out.println("p >>>>> "+" "+l);
//            for(int i=0;i<n;i++)
//                System.out.println(p[i]);
                
        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p[i]=vv[i].no;
        	performanceTask[l]=p;
    	}
    	
    	
    	
//   ------------------------------------------------------------------------
        for(int l=0;l<noOfFN;l++){
       		int n=noOfTask;
        	int c=2;
            double[][] B=new double[n][2];
            B=decisionMatrixFN[l];

        	//Normalize the decision matrix using Eq. (19).
        	for(int j=0;j<c;j++){
            	double square_sum=0;
            	for(int i=0;i<n;i++) square_sum+=B[i][j]*B[i][j];
        	    square_sum=Math.sqrt(square_sum);
       	 		for(int i=0;i<n;i++) B[i][j]=(B[i][j]/square_sum); // Eq. (19).
        	}

        	// Calculate the weighted normalized decision matrix based on Eq. (20)
        	for(int j=0;j<c;j++){
            	for(int i=0;i<n;i++) B[i][j]=B[i][j]*weightMatrixFN[l][j];  // Eq. (20)
        	}


        	double[] positive=new double[c];
        	double[] negative=new double[c];
        	// Evaluate the positive and negative ideal solutions for each criterion k as 
        	// minimum and maximum value of kth column weighted normalized matrix.
        	for(int j=0;j<c;j++){
            	negative[j]=positive[j]=B[0][j];
            	for(int i=0;i<n;i++){
                	negative[j]=max(negative[j],B[i][j]);
                	positive[j]=min(positive[j],B[i][j]);
            	}
        	}

        	double[] d_pos=new double[n];
        	double[] d_neg=new double[n];
        	double[] p=new double[n];
 			// Determine the distance of each alternative from the positive and negative ideal 
 			// solutions based on Eq.(21) and (22).
        	for(int i=0;i<n;i++){
            	for(int j=0;j<c;j++){
                	d_pos[i]+=(B[i][j]-positive[j])*(B[i][j]-positive[j]); 
                	d_neg[i]+=(B[i][j]-negative[j])*(B[i][j]-negative[j]);   
            	}
            	d_pos[i]=Math.sqrt(d_pos[i]);  // Eq.(21)
            	d_neg[i]=Math.sqrt(d_neg[i]);  // Eq.(22)
                                 
        	}

        	// Compute the performance score for each alternative following Eq. (23),
        	for(int i=0;i<n;i++){
        	 	p[i]=d_neg[i]/(d_neg[i]+d_pos[i]);     // Eq.(23)    
        	}

        	// ranked in decreasing order of their performance scores 
        	//Vector<Pair <float, Integer>> vv=new Vector<Pair <float, Integer>>(n);
        	PAIR[] vv=new PAIR[n];
        	for(int i=0;i<n;i++){
        		vv[i]=new PAIR(p[i],i);
        	}
        	//sort_reverse(vv,n);
        	for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    PAIR temp = new PAIR(0.0,0);
                    if (vv[j].val < vv[i].val) {
                        // Swapping
                        temp = vv[i];
                        vv[i] = vv[j];
                        vv[j] = temp;
                    }
                }
            }
            PAIR temp = new PAIR(0.0,0);
        	int i1=0,j1=n-1;
        	while(i1<=j1){
        	    temp = vv[i1];
                vv[i1] = vv[j1];
                vv[j1] = temp;
                i1++;
                j1--;
        	}
        	for(int i=0;i<n;i++)
        		p[i]=vv[i].no;
        	performanceFN[l]=p;
    	}
    
    }


//    public static void matching(){
//	    // Q -- Quota of each FN
//	   	Integer[] Q=new Integer[noOfFN];
//	   	for(int i=0;i<noOfFN;i++)
//	   	    Q[i]=50; // store some value
//	
//	    // Assign  -- Assigned tasks in FN
//	    Integer[][] Assign=new Integer[noOfFN][noOfTask];
//	    	for(int i=0;i<noOfTask;i++){
//	    	    double[] tj=new double[noOfFN];
//	    	    tj=performanceTask[i];
//	        	int n=noOfFN;
//	        	for(int x=0;x<n;x++){
//	        		// fi∗= highest ranked FN in P(tj) to which tj has not proposed yet
//	        		// Send proposal to fi∗.
//	            	int fi=(int)tj[x];
//	            	if(Q[fi]>0){    // if Qi∗ > 0 then
//	                	Assign[fi][x]=i+100;
//	//                	System.out.println("Task "+i+" is being assigned to Fog Device "+fi);
//	                	Q[fi]=Q[fi]-1;
//	                	break;
//	            	}
//	            	else{
//	                	// Reject the assignment request;
//	            		count++;
//	           	}
//	        }	    
//	    }
//	}

	public static void matching(){
			
	    // Assign  -- Assigned tasks in FN
	    Integer[][] Assign=new Integer[noOfFN][noOfTask];
	    int count = 0; // Initialize the outage count to zero
	    
	    for(int i=0;i<noOfTask;i++){
	        double[] tj=new double[noOfFN];
	        tj=performanceTask[i];
	        int n=noOfFN;
	        int devices = numberOfDependentDevices[i];
	        boolean assigned = false; // Initialize the assigned flag to false
	        for(int x=0;x<n;x++){
	            // fi∗= highest ranked FN in P(tj) to which tj has not proposed yet
	            // Send proposal to fi∗.
	            int fi=(int)tj[x];
	            if(Q[fi]>0){    // if Qi∗ > 0 then
	                if(Q[fi] > devices ) {			//required VRU is less than Quota of current task
	                	Assign[fi][x]=i+100;
	                	Q[fi]=Q[fi]-devices;
	                	assigned = true; // Set the assigned flag to true
	                	break;
	                }
	            }
	        }
	        if(!assigned){ // If the task was not assigned, increment the outage count
	            count++;
	        }
	    }
	    System.out.println("Total number of outages: " + count); // Print the total number of outages
	}


}
