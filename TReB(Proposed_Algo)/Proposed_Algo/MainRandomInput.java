package Proposed_Algo;

import java.util.HashMap;
import java.util.Scanner;

import Proposed_Algo.Main.Binarycases;
import Proposed_Algo.Main.criteria_Weight;
import Proposed_Algo.Main.fogNode;
import Proposed_Algo.Main.iotDevice;

public class MainRandomInput {
	static int caseDistribution[] = new int[4];
	
	static class criteria_Weight
	  {
	    double CW_of_Energy_of_iot;
	    double CW_of_Latency_of_iot;
	    double CW_of_Computational_demand_of_iot;
	    double CW_of_Energy_of_fog;
	    double CW_of_deadline_of_fog;
	  }
	  //class for fog node ranking using AHP
	  static class fogNode
	  {
	    int Id;
	    double computational_power_of_fogNode;
	    double power_of_fogNode;
	    double deadline;
	    double computational_demand;
	    double energy;
	    double inputSize;
	    double weightage;
	    double quota;

	  }
	  static class iotDevice
	  {
	    int Id;
	    double computational_demand;
	    double inputSize;
	    double energy;
	    double latency;
	    double weightage;
	    double dependency;
	  }
	  static class Binarycases
	  {
	    double energy;
	    double total_time;
	  }

	  public static double log2 (double x)
	  {
	    return Math.log (x) / Math.log (2);
	  }
	  public static int MatchingAlgo (iotDevice[]iot, fogNode[]fog, int count)
	  {
	    //we are implementing one to many like every iot task will be having multiple options to go 
	    int total_number_of_iot_tasks = iot.length;
	    //storing the matching of iot, fog 
	    int total_number_of_fog_nodes = fog.length;
	    HashMap < Integer, Integer > matching = new HashMap < Integer, Integer > ();
	    //to check weatrher the task went to any fog device or any fog assigned or not to that particular task
	    int flag = 0;
	    for (int i = 0; i < total_number_of_iot_tasks; i++)
		{
			flag = 0;		//reinitializing the flag 0 to track again
		
			for (int j = 0; j < total_number_of_fog_nodes; j++)
			  {
			    if (fog[j].quota >= iot[i].dependency)
			      {
				fog[j].quota = fog[j].quota - iot[i].dependency;
				matching.put (iot[i].Id, fog[j].Id);
				flag = 1;	//to track that is the iot task assigned to a fog node or node
		
		
			      }
			    if (flag == 1)	//if task is assigned according to preference then breaking the loop and going on next task
			      break;
		
			  }
			if (flag == 0)
			  {
			    matching.put (iot[i].Id, -1);
			    count++;
			  }
	      }





	    //printing pairs of HashMap
//	    System.out.println ("printing matched iot fog pair ...");
	//  for (Map.Entry m:matching.entrySet ())
//	      {
//		if ((Integer) m.getValue () != -1)	//if none of the fog node is assigned then  it will go in the else case 
//		  System.out.println ("iot task number " + m.getKey () +
//				      " is matched with " +
//				      "fog node " + m.getValue ());
//		else
//		  System.out.println ("iot task number" + m.getKey () +
//				      " is not matched with any fog node due to busyness of fog nodes");
//	      }

	    return count;
	  }
	  
	  
	  public static criteria_Weight AHP ()
	  {
	    criteria_Weight cw = new criteria_Weight ();

	    //for fog node 
	    //oth row/column energy
	    //1st row/column latency
	    //2nd  row/column computational demand
	    //for fog node 
	    double sum = 0;
	    int n = 3;
	    //for fog node
	    double[][] standard_matrix = { {1.0, 4.0, 5.0},
	    {0.25, 1.0, 3.0},
	    {0.2, 0.33, 1.0}
	    };
	    //step one to calculate C.W.
	    double[][] matrix_for_CW = new double[n][n];
	    for (int i = 0; i < n; i++)
	      {
		sum = 0;
		for (int j = 0; j < n; j++)
		  {
		    sum = sum + standard_matrix[j][i];

		  }
		for (int j = 0; j < n; j++)
		  {
		    matrix_for_CW[j][i] = standard_matrix[j][i] / sum;
		  }

	      }

	    //since first is energy then latency then computational demand
	    sum = 0;

	    for (int i = 0; i < n; i++)
	      {
		sum = sum + matrix_for_CW[0][i];
	      }

	    cw.CW_of_Energy_of_iot = sum / n;
	    sum = 0;
	    for (int i = 0; i < n; i++)
	      {
		sum = sum + matrix_for_CW[1][i];
	      }
	    cw.CW_of_Latency_of_iot = sum / n;
	    sum = 0;
	    for (int i = 0; i < n; i++)
	      {
		sum = sum + matrix_for_CW[2][i];
	      }
	    cw.CW_of_Computational_demand_of_iot = sum / n;



	    //matrix for WS 
	    double[][] matrix_for_WS = new double[n][n];
	    //defining variables for weighted sum accordingly
	    double weighted_sum_of_energy = 0;
	    double weighted_sum_of_latency = 0;
	    double weighted_sum_of_computational_demand = 0;
	    sum = 0;
	    for (int i = 0; i < n; i++)
	      {
		matrix_for_WS[i][0] = cw.CW_of_Energy_of_iot * standard_matrix[i][0];

	      }
	    for (int i = 0; i < n; i++)
	      {
		matrix_for_WS[i][1] = cw.CW_of_Latency_of_iot * standard_matrix[i][1];

	      }
	    for (int i = 0; i < n; i++)
	      {
		matrix_for_WS[i][2] =
		  cw.CW_of_Computational_demand_of_iot * standard_matrix[i][2];

	      }

	    //calculating weighted sum for all
	    sum = 0;

	    for (int i = 0; i < n; i++)
	      {
		weighted_sum_of_energy = weighted_sum_of_energy + matrix_for_CW[0][i];
	      }

	    for (int i = 0; i < n; i++)
	      {
		weighted_sum_of_latency =
		  weighted_sum_of_latency + matrix_for_CW[1][i];
	      }

	    for (int i = 0; i < n; i++)
	      {
		weighted_sum_of_computational_demand =
		  weighted_sum_of_computational_demand + matrix_for_CW[2][i];
	      }

	    double Energy_ratio = weighted_sum_of_energy / cw.CW_of_Energy_of_iot;

	    double Latency_ratio = weighted_sum_of_latency / cw.CW_of_Latency_of_iot;
	    double Computational_demand_ratio =
	      weighted_sum_of_computational_demand /
	      cw.CW_of_Computational_demand_of_iot;
	    double lambda_max =
	      (Energy_ratio + Latency_ratio + Computational_demand_ratio) / 3;
	    double CI = (lambda_max - n) / (n - 1);
	    double CR = CI / 0.58;	//random index is 0.58 for n=3;
	    if (CR < 0.10)
	      {
		System.out.println ("Matrix is consistent");
	      }

	    //cw for iot device
	    //oth row/column energy
	    //1st row/column deadline
	    double[][] standard_matrix_for_iot = { {1.0, 5.0},
	    {0.2, 1.0}
	    };
	    n = 2;
	    double[][] matrix_for_CW_of_iot = new double[n][n];
	    for (int i = 0; i < n; i++)
	      {
		sum = 0;
		for (int j = 0; j < n; j++)
		  {
		    sum = sum + standard_matrix_for_iot[j][i];

		  }
		for (int j = 0; j < n; j++)
		  {
		    matrix_for_CW_of_iot[j][i] = standard_matrix_for_iot[j][i] / sum;
		  }

	      }
	    sum = 0;
	    for (int i = 0; i < n; i++)
	      {
		sum = sum + matrix_for_CW_of_iot[0][i];
	      }

	    cw.CW_of_Energy_of_fog = sum / n;
	    sum = 0;
	    for (int i = 0; i < n; i++)
	      {
		sum = sum + matrix_for_CW_of_iot[1][i];
	      }
	    cw.CW_of_deadline_of_fog = sum / n;

	    return cw;


	  }


	  public static fogNode[] fogpreference (double downlink_data_rate,
						 double uplink_data_rate,
						 double output_size,
						 double capacity_of_VRU)
	  {
	    fogNode[]fog;
	    int number_of_fogNode;
	    Scanner sc = new Scanner (System.in);	//System.in is a standard input stream  
	    System.out.println ("Enter the number of fogNodes");
	    int n = sc.nextInt ();
	    fog = new fogNode[n];
	    for (int i = 0; i < n; i++){
			fog[i] = new fogNode ();
			fog[i].Id = i + 1;
			System.out.println ("enter the following for fog node:" + (i + 1));
			System.out.println ("power of fognode");
			fog[i].power_of_fogNode = sc.nextDouble ();
			System.out.println ("computational power  of fognode");
			fog[i].computational_power_of_fogNode = sc.nextDouble ();
			System.out.println ("deadline of fognode");
			fog[i].deadline = sc.nextDouble ();
			System.out.println ("inputSize of fognode in Kb");
			fog[i].inputSize = sc.nextDouble ();
			System.out.println ("computational_demand of fognode");
			fog[i].computational_demand = sc.nextDouble ();
			System.out.println ("enter quota of the fognde:: ");
			fog[i].quota = sc.nextDouble ();
		
			double upload_latency = fog[i].inputSize / (uplink_data_rate * 1000);
			double download_latency = output_size / (downlink_data_rate * 1000);
			double computational_latency = (fog[i].computational_demand ) / capacity_of_VRU;
			criteria_Weight cw = AHP ();
			fog[i].energy = fog[i].power_of_fogNode * (upload_latency + download_latency) + computational_latency * fog[i].computational_power_of_fogNode;
			fog[i].weightage = cw.CW_of_Energy_of_fog * fog[i].energy +  cw.CW_of_deadline_of_fog * fog[i].deadline;
	    }
	    // sorting the fog nodes based on preference 
	    for (int i = 0; i < n; i++){
		for (int j = i + 1; j < n; j++){
		    if (fog[i].weightage < fog[j].weightage)
		      {
				fogNode temp = fog[i];
				fog[i] = fog[j];
				fog[j] = temp;
		      }
		  }
	    }
	    System.out.println ("order of the fog nodes are following\n");
	    for (int i = 0; i < n; i++)
	      {
	    	System.out.print (fog[i].Id + " ");
	      }
	    return fog;
	  }

	  //function for both local execution
	  public static Binarycases both_locally_executed (double
							   local_execution_time,
							   double uplink_time,
							   double downlink_time,
							   double power_of_iotd1,
							   double power_of_iotd2,
							   double
							   computational_power_of_fognode,
							   double computational_delay,
							   double X, double S,
							   double N,
							   double power_of_fogDevice)
	  {
	    double total_time =
	      2 * local_execution_time + uplink_time + downlink_time;
	    double total_energy =
	      power_of_iotd1 * (local_execution_time) +
	      power_of_iotd2 * (local_execution_time) + (uplink_time +
							 downlink_time) *
	      power_of_fogDevice;

	    double computation_time_for_iotd1 =
	      (uplink_time + computational_delay) * (X - 1);
	    double computation_time_for_iotd2 =
	      (uplink_time + computational_delay) * (S - 1);
	    double remaining_time = (uplink_time + computational_delay) * (N - S);
	    double waiting_time =
	      Math.max (computation_time_for_iotd1, computation_time_for_iotd2);
	    total_energy =
	      total_energy + (X - 2 +
			      N) * (computational_power_of_fognode *
				    computational_delay +
				    power_of_fogDevice * (uplink_time +
							  downlink_time));

	    total_time = total_time + remaining_time + waiting_time;
	    Binarycases bn = new Binarycases ();
	    bn.total_time = total_time;
	    bn.energy = total_energy;

	    System.out.println ("CASE 1: both tasks are executed locally!!");
	    // System.out.println("total time is:: "+total_time);
	    //     System.out.println("total energy is:: "+total_energy+"\n");
	    return bn;

	  }

	  //function for device 1 on fog node and 2 local execution
	  public static Binarycases first_onfog_second_locally (double
								local_execution_time,
								double uplink_time,
								double downlink_time,
								double power_of_iotd1,
								double power_of_iotd2,
								double
								computational_power_of_fognode,
								double
								computational_delay,
								double X, double S,
								double N,
								double
								power_of_fogDevice)
	  {
	    double total_time =
	      local_execution_time + uplink_time + downlink_time +
	      computational_delay;
	    double total_energy =
	      power_of_fogDevice * (uplink_time + downlink_time) +
	      computational_power_of_fognode * computational_delay +
	      power_of_iotd2 * (local_execution_time);
	    double computation_time_for_iotd1 =
	      (uplink_time + computational_delay) * (X - 1);
	    double computation_time_for_iotd2 =
	      (uplink_time + computational_delay) * (S - 1);
	    double remaining_time = (uplink_time + computational_delay) * (N - S);
	    double waiting_time =
	      Math.max (computation_time_for_iotd1, computation_time_for_iotd2);
	    total_energy =
	      total_energy + (X - 2 +
			      N) * (computational_power_of_fognode *
				    computational_delay +
				    power_of_fogDevice * (uplink_time +
							  downlink_time));


	    total_time = total_time + remaining_time + waiting_time;
	    Binarycases bn = new Binarycases ();
	    bn.total_time = total_time;
	    bn.energy = total_energy;
	    System.out.println ("CASE 2: iotd1 execute at fog and d2 locally!!");
	    // System.out.println("total time is:: "+total_time);
	    //     System.out.println("total energy is:: "+total_energy+"\n");
	    return bn;

	  }

	  //function for task 1 local execution(on iotd1) and  task 2 on fog node
	  public static Binarycases first_locally_second_onfog (double
								local_execution_time,
								double uplink_time,
								double downlink_time,
								double power_of_iotd1,
								double power_of_iotd2,
								double
								computational_power_of_fognode,
								double
								computational_delay,
								double X, double S,
								double N,
								double
								power_of_fogDevice)
	  {
	    double total_time =
	      local_execution_time + 2 * uplink_time + computational_delay;
	    double total_energy =
	      (2 * uplink_time * power_of_fogDevice) +
	      computational_delay * computational_power_of_fognode +
	      power_of_iotd1 * (local_execution_time);
	    double computation_time_for_iotd1 =
	      (uplink_time + computational_delay) * (X - 1);
	    double computation_time_for_iotd2 =
	      (uplink_time + computational_delay) * (S - 1);
	    double remaining_time = (uplink_time + computational_delay) * (N - S);
	    double waiting_time =
	      Math.max (computation_time_for_iotd1, computation_time_for_iotd2);
	    total_energy = total_energy + (X - 2 + N) * (computational_power_of_fognode * computational_delay + power_of_fogDevice * (uplink_time + downlink_time));	//X and Sth task are dependent so X-1+S-1+Y-S will be remaining task that is X-2+N

	    total_time = total_time + remaining_time + waiting_time;
	    Binarycases bn = new Binarycases ();
	    bn.total_time = total_time;
	    bn.energy = total_energy;
	    System.out.println ("CASE 3: iotd1 execute locally and d2 at fog!!");
	    // System.out.println("total time is:: "+total_time);
	    //     System.out.println("total energy is:: "+total_energy+"\n");

	    return bn;
	  }
	  //both execute on fog node (iotd1 and iotd2)
	  public static Binarycases both_onfog (double local_execution_time,
						double uplink_time,
						double downlink_time,
						double power_of_iotd1,
						double power_of_iotd2,
						double computational_power_of_fognode,
						double computational_delay, double X,
						double S, double N,
						double power_of_fogDevice)
	  {
	    double total_time = 2 * uplink_time + 2 * computational_delay;
	    double total_energy =
	      computational_power_of_fognode * (2 * uplink_time +
						2 * computational_delay);
	    double computation_time_for_iotd1 =
	      (uplink_time + computational_delay) * (X - 1);
	    double computation_time_for_iotd2 =
	      (uplink_time + computational_delay) * (S - 1);
	    double remaining_time = (uplink_time + computational_delay) * (N - S);
	    double waiting_time =
	      Math.max (computation_time_for_iotd1, computation_time_for_iotd2);
	    total_energy =
	      total_energy + (X - 2 +
			      N) * (computational_power_of_fognode *
				    computational_delay +
				    power_of_fogDevice * (uplink_time +
							  downlink_time));


	    total_time = total_time + remaining_time + waiting_time;
	    Binarycases bn = new Binarycases ();
	    bn.total_time = total_time;
	    bn.energy = total_energy;
	    System.out.println ("CASE 4: iotd1 and d2 both execute at fog!!");
	    // System.out.println("total time is:: "+total_time);
	    //     System.out.println("total energy is:: "+total_energy+"\n");

	    return bn;
	  }


	  public static iotDevice[] iotpreference (double local_execution_time,
						   double uplink_time,
						   double downlink_time,
						   double power_of_iotd1,
						   double power_of_iotd2,
						   double computational_power_of_fognode,
						   double computational_delay,
						   double X, double S, double N,
						   double power_of_fogDevice,
						   double downlink_data_rate,
						   double uplink_data_rate,
						   double output_size,
						   double capacity_of_VRU)
	  {
	    Scanner sc = new Scanner (System.in);
	    iotDevice[]iot;
	    
	    System.out.println ("\nEnter the number of independent task : ");
	    int number_of_indepenenttask = sc.nextInt ();
	    //4 is used for binary cases 
	    int n = number_of_indepenenttask + 4;
	    iot = new iotDevice[n];
	    Binarycases bn = new Binarycases ();

	    criteria_Weight cw = AHP ();
	    
	    int range = 4 - 1 + 1;
	    
	    for (int i = 0; i < 4; i++)
	      {
	    	int random = (int)(Math.random() * range);
	    	caseDistribution[i] = random+1;
	    	
			if (random == 0)
			  {
			    bn = both_locally_executed (local_execution_time, uplink_time,
						     downlink_time, power_of_iotd1,
						     power_of_iotd2,
						     computational_power_of_fognode,
						     computational_delay, X, S, N,
						     power_of_fogDevice);
		
			  }
			else if (random == 1)
			  {
			    // else if(iotd1_execute_onfognode==1&&iotd2_execute_onfognode==0)
			    bn = first_onfog_second_locally (local_execution_time, uplink_time,
							  downlink_time, power_of_iotd1,
							  power_of_iotd2,
							  computational_power_of_fognode,
							  computational_delay, X, S, N,
							  power_of_fogDevice);
		
			  }
			else if (random == 2)
			  {
			    // else if(iotd1_execute_onfognode==0&&iotd2_execute_onfognode==1)
			    bn = first_locally_second_onfog (local_execution_time, uplink_time,
							  downlink_time, power_of_iotd1,
							  power_of_iotd2,
							  computational_power_of_fognode,
							  computational_delay, X, S, N,
							  power_of_fogDevice);
		
			  }
			else
			  {
			    // else 
			    bn = both_onfog (local_execution_time, uplink_time, downlink_time,
					  power_of_iotd1, power_of_iotd2,
					  computational_power_of_fognode, computational_delay,
					  X, S, N, power_of_fogDevice);
			  }
	   		  iot[i] = new iotDevice ();
			  iot[i].Id = i + 1;
			  iot[i].energy = bn.energy;
			  iot[i].dependency = 2;
			  iot[i].latency = bn.total_time;
			  iot[i].computational_demand = 210;
			  iot[i].inputSize = 400;
			  iot[i].weightage =  cw.CW_of_Energy_of_iot * iot[i].energy +
					  	cw.CW_of_Computational_demand_of_iot * iot[i].computational_demand +cw.CW_of_Latency_of_iot * iot[i].latency;
	      }
	    
	    
	    for (int i = 4; i < n; i++)
		{
			iot[i] = new iotDevice ();
			iot[i].Id = i + 1;
			iot[i].dependency = 1;
			System.out.println ("enter computational demand of  task :: " + (i + 1));
			iot[i].computational_demand = sc.nextDouble ();
			System.out.println ("enter inputSize of " + " task:: " + (i + 1));
			iot[i].inputSize = sc.nextDouble ();
			if (iot[i].computational_demand < 300)	//if computational demand is less than 300 then it'll be executed locally
			  {
			    iot[i].latency = iot[i].computational_demand / 16;
			    iot[i].energy = power_of_iotd1 * iot[i].latency;
			  }
			else
			  {
			    //executing on fog node
			    iot[i].latency = (iot[i].inputSize / (uplink_data_rate * 1000) +	//upload
					      (iot[i].computational_demand) /  capacity_of_VRU +	//computational
					      (output_size / (downlink_data_rate * 1000)));	//download
			    iot[i].energy = power_of_fogDevice * (iot[i].inputSize /
						    (uplink_data_rate * 1000) +
						    output_size / (downlink_data_rate * 1000)) +
			      computational_power_of_fognode * (iot[i].computational_demand) /  capacity_of_VRU;
			    
			  }
			//calculating weightage for ranking purpose
			iot[i].weightage =
			  cw.CW_of_Energy_of_iot * iot[i].energy +
			  cw.CW_of_Computational_demand_of_iot * iot[i].computational_demand +
			  cw.CW_of_Latency_of_iot * iot[i].latency;
		 }

	    //ordering the iot devices according to their weighttage of criteria weights
	    for (int i = 0; i < n; i++)
		  {
			for (int j = i + 1; j < n; j++)
			  {
			    if (iot[i].weightage < iot[j].weightage)
			      {
				iotDevice temp = iot[i];
				iot[i] = iot[j];
				iot[j] = temp;
			      }
			  }
		   }
	    
	    System.out.println ("order of the iot devices are following\n");
	    for (int i = 0; i < n; i++)
	      {
	    	System.out.print (iot[i].Id + " ");
	      }
	    
	    System.out.println ("\n");
	    return iot;
	  }


	  public static double log10 (int x)
	  {
	    return Math.log (x) / Math.log (10);
	  }

	  public static void main (String[]args)
	  {
	    Scanner sc = new Scanner (System.in);
	    //taking all required inputs
	    System.out.println ("please provide the following inputs::\n");
	    System.out.println ("enter the computational demand of the task in million cycle::");
	    double computational_demand = 300; //sc.nextDouble();
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
	    double output_size = 15; //sc.nextDouble ();
	    System.out.println ("enter the computational power of fog node::");
	    double computational_power_of_fognode = 0.35; //sc.nextDouble ();

	    System.out.println ("enter the value of X::");
	    double X = 3; //sc.nextDouble ();
	    System.out.println ("enter the value of S::");
	    double S = 4; //sc.nextDouble ();
	    System.out.println ("enter the value of N::");
	    double N = 5; //sc.nextDouble ();


	    double local_execution_time = computational_demand / frequency;
	    System.out.println ("local execution time is:: " + local_execution_time + "\n");

	    double log10distance_to_fog_node = log10 (distance_to_fog_node);

	    double PLIJ = 38.02 + 20 * log10distance_to_fog_node;
	    //   System.out.println(PLIJ);
	    double power_to_calculate = -(PLIJ / 10);
	    double channel_gain = Math.pow (10, power_to_calculate);

	    double noise_power = Math.pow (10, -10);

	    double in_the_log = 1 + (power_of_iot_device * channel_gain) / noise_power;

	    double log2in_the_log = log2 (in_the_log);
	    double uplink_data_rate = bandwidth * log2in_the_log;

	    System.out.println ("uplink data rate is in megahetz :: " + uplink_data_rate + "\n");

	    double downlink_data_rate = bandwidth * log2 (1 + (power_of_fogDevice * channel_gain) / noise_power);
	    System.out.println ("downlink  data rate is in megahetz :: " + downlink_data_rate + "\n");
	    double Energy_of_IOT_device = power_of_iot_device * local_execution_time;
	    System.out.println ("enery of iot device is in joule :: " +	Energy_of_IOT_device + "\n");
	    
//	    double computational_capacity_of_one_fognode = 50000;
//	    double Iij = 300;
	// 
//	    double computational_delay =
//	      (computational_demand * Iij) / (computational_capacity_of_one_fognode);
	    
	    double fog_node_capacity = 6 * Math.pow(10, 9);
	    double temp_quota = 50.0;
	    double capacity_of_VRU = fog_node_capacity/temp_quota ;
	    
	    double computational_delay = (computational_demand) / (capacity_of_VRU);

	    double uplink_data_rate_time = output_size / (uplink_data_rate * 1000);
	    double downlink_data_rate_time = output_size / (downlink_data_rate * 1000);



	    //calling function to rank fognodes 
	    fogNode[]fog;
	    iotDevice[]iot;
	    //calling fog preferrence to get the order of fog nodes 
	    fog = fogpreference (downlink_data_rate, uplink_data_rate, output_size, capacity_of_VRU);
	    //calling iot preferrence to get the order of iot devices
	    iot = iotpreference (local_execution_time, uplink_data_rate_time,
			     downlink_data_rate_time, power_of_iot_device,
			     power_of_iot_device, computational_power_of_fognode,
			     computational_delay, X, S, N, power_of_fogDevice,
			     downlink_data_rate, uplink_data_rate, output_size, capacity_of_VRU);
	    
	    //Calling Matching Theory algorithm
	    
	    int count = MatchingAlgo (iot, fog, 0);


	    //printing total energies
//	    System.out.println("Printing Energy values for Fog Nodes  ---------->");
	    double total_energy = 0.0;
	    int no_fog_nodes = fog.length;
	    for(int i=0; i<no_fog_nodes; i++) {
	    	total_energy += fog[i].energy;
	    }

	    int no_iots = iot.length;
	    System.out.println("No of IoT : "+no_iots);
	    for(int i=0; i<no_iots; i++) {
	    	total_energy += iot[i].energy;
	    }
	    System.out.println("\nTotal Energy for " + no_iots +" devices = "+ total_energy);
	    
	    double total_latency = 0.0;
	    for(int i=0; i<no_iots; i++) {
	    	total_latency += iot[i].latency;
	    }
	    System.out.println("\nTotal Latency for " + no_iots +" devices = "+ total_latency);
	    
	    System.out.println("Number of Outages : " + count);
	    
	    for(int i=0; i<4; i++) {
	    	System.out.print( caseDistribution[i] + " ");
	    }
	    
	    System.out.println();
	    sc.close();
	  }
	 
}
