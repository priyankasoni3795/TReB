# TReB
# Execution Environment:
Operation System: Ubuntu 22.04.1, 64bit. 

Physical Memory (RAM) 32 GB.

# Prerequisites
 Java SE Development Kit 19.0.2
 
 Eclipse IDE 64 bit for Java Developers 

# Installation
 -> Download Eclipse IDE 64 bit from https://www.eclipse.org/downloads/packages/release/kepler/sr1/eclipse-ide-java-developers.
 
 -> Then extract the downloded zip folder and Open eclipse-installer in terminal.
 
 -> Then install Eclipse by typing   ./eclipse-inst
 
 -> Now install the Eclipse IDE for java developer, after that Eclipse IDE will visible and it automatically took the JAVA version and you can browse for latest version of JAVA 19.0.2 
 
 -> It will install Eclipse under the following location and then click install.
 
 -> Now set the license agreement and click accept now, it will successfully instal Eclipse.
 
 -> Now click launch to open it.


# Code Execution 
 -> There is Main.java file which contains the code of proposed algorithm TReB for uniform generation. We need to to simply run that java file to get the output. 

 -> There is a MatchingAlgo function that implements a one-to-many matching theory algorithm for assigning each IoT task to a single fog node.

 -> There is a criteria_Weight AHP function that implements the Analytic Hierarchy Process algorithm by evaluating criteria weights for task preference.
 
 -> Once execution started, it will ask for some properties of components used in our architecture such as double computational_demand U[210-480] Mcycles, IoT device frequency (16)MHz, channel frequency(bandwidth) (10)MHz , Distance of IoT device to fog node (100)m, power of iot device U[0.1 to 1]W, power of fog node U[1 to 2]W, output size U[10,20]Kb, computational power of fog node U[0.35, 0.55]W. These values we have defined the range and to make execution simple we have assigned the pre-determined values to these components. 
 
 -> So, after expecting the code, it will take that pre-determined values from code, then it will ask for the details related to fogNodes. The details such as computational power of fogNode, power of fogNode, deadline, computational demand, energy, inputSize, weightage, quota.
 
-> Once consumed the details of fog nodes, the program will ask for details related to tasks. It will first ask for number of independent tasks. Then for each task it will ask to enter computational demand of each task and input size of each subtask. 

-> For simplicity we have given a text file which contains sample input values for number of independent tasks, computational demand of and input size of each subtask. First line this text file represents number of independent task. Then from second line, each pair of two lines represent computational demand and input size. These pair of values repeat for given number of tasks.

-> The sample input files are provided for both homogeneous environment and heterogeneous environment.

-> Now all the inputs are feeded to the program, it will perform the calculation as per proposed algorithm TReB and will generate the output. The output values given by program are Total Energy, Total Latency and Number of Outages(if any).

-> There is a Dependent.java file which contains the code of proposed algorithm TReB for random generation. We need to to simply run that java file to get the output.

# Usage

-> For both Uniform and Random generation 

-> We have taken the number of IoT devices that vary in the range of [250, 500, 750, 1000] for both homogeneous and heterogeneous scenarios.

-> We have set 5 Fog nodes in the fog netwrok for both homogeneous and heterogeneous scenarios. 

-> We carry out 10 test runs for every scenario. We have added the input file for each parameter.

-> Important Parameters such as - 

-> Input size of task is 300 Kb U[300-600] Kb

-> Output size of task is 10 Kb U[10-20] Kb

-> Arduino frequency is 16MHz for local task execution.

-> Bandwidth is 10MHz 

-> Computational demand of task is U[210,480] Mcycles which will decide the task execution.

-> Maximum transmission power of an IoT device is  U[0.1 to 1] W

-> Maximum transmission power of a FN  U[1 to 2] W

-> Computational power of FN  U[0.35,0.55] W

-> Computation rate of FN  U[6,10] GHz

-> Deadline of a task  U[30-60] s

-> Quota of FN  [50-500]

# References
[1] C. Swain, M. N. Sahoo, A. Satpathy, K. Muhammad, S. Bakshi, J. J. Rodrigues, V. H. C. de Albuquerque, METO:
Matching-theory-based efficient task offloading in IoT-fog interconnection networks, IEEE Internet of Things Jour-
nal 8 (16) (2020) 12705–12715. doi: https://doi.org/10.1109/JIOT.2020.3025631

[2] J. Yan, S. Bi, Y. J. Zhang, M. Tao, Optimal task offloading and resource allocation in mobile-edge computing
with inter-user task dependency, IEEE Transactions on Wireless Communications 19 (1) (2019) 235–250. doi: https://doi.org/10.1109/TWC.2019.2943563

[3] F. Chiti, R. Fantacci, B. Picano, A matching theory framework for tasks offloading in fog computing for IoT
systems, IEEE Internet of Things Journal 5 (6) (2018) 5089–5096. doi: https://doi.org/10.1109/JIOT.2018.2871251

# Contributors

-> Miss Priyanka Soni

   https://scholar.google.com/citations?user=LZKL3o4AAAAJ&hl=en

-> Mr Ajay Gajanan Hajare

  https://github.com/AjayHajare

-> Mr. Keerthan Kumar T G

   https://scholar.google.com/citations?user=fW7bzK8AAAAJ&hl=en

-> Dr. Sourav kanti addya

   https://souravkaddya.in/

# Contact

If you have any questions, simply write a mail to sonipriyanka31994(AT)gmail(DOT)com



