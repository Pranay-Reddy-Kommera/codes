#include "TS_Data_GPU.h"

int main()
{
	Performance pm;

	//Variables
	string *angleFileName = new string[NUM_BUSES]; //string object to store the file name
	string *powerFileName = new string[NUM_BUSES]; //string object to store the file name
	float *PMU_Angle = new float[TOTAL_PMU*SAMPLES];
	float *PMU_Power = new float[TOTAL_PMU*SAMPLES];
	float *angle_Samples = new float[TOTAL_PMU];
	float *power_Samples = new float[TOTAL_PMU];
	float *centerOfAngle = new float[1];
	float *product_PMU = new float[TOTAL_PMU];

	pm.OnSetStartTime();
	//Function to generate file names
	GenerateFileNames(angleFileName, powerFileName);

	//Function to print file names
	//If you want to check if file names are assigned properly uncomment the line below
	//PrintFileNames(angleFileName, powerFileName);

	//Function to extract data from the files
	ExtractData(angleFileName, powerFileName, PMU_Angle, PMU_Power);

	//cout << "Angle" << PMU_Angle[362] << endl;
	//cout << "Angle" << PMU_Power[362] << endl;

	PMU_Angle[1089] = 1200;
	PMU_Angle[1087] = 1200;


	//Extracting data for COA
	for (int p = 0; p < TOTAL_PMU; p++)
	{
		angle_Samples[p] = PMU_Angle[p];
		power_Samples[p] = PMU_Power[p];
	}

	pm.OnSetStopTime();
	float time_Data_Extraction = pm.OnGetElapsedTime();
	
	
	//GPU Function to perform center of angles
	float time_COA=CenterOfAnglesGPU(pm, angle_Samples, power_Samples,centerOfAngle, product_PMU);
	

	//CPU Function to extract center of angle
	pm.OnSetStartTime();
	CenterOfAngle(angle_Samples, power_Samples,centerOfAngle);
	pm.OnSetStopTime();
	float time_COA_CPU = pm.OnGetElapsedTime();

	

	//cout << "The center of angle is " << centerOfAngle[0] << endl;

	//GPU Function to find stability 
	float time_Stability=StabilityGPU(pm,PMU_Angle, centerOfAngle);
	
	float total_Time_Without_Extraction_CPU = time_COA_CPU + time_Stability;
	float total_Time_With_Extraction_CPU = (time_Data_Extraction*1000) + time_COA_CPU + time_Stability;
	float total_Time_Without_Extraction_GPU = time_COA + time_Stability;
	float total_Time_With_Extraction_GPU = (time_Data_Extraction*1000) + time_COA + time_Stability;

	cout << "Data extraction time " << time_Data_Extraction * 1000 << endl;
	cout << "The total time for GPU COA is " << time_COA << endl;
	cout << "The time taken for CPU_COA is " << time_COA_CPU * 1000 << endl;
	cout << "The time taken for CPU_Stability is " << time_Stability << endl;

	cout << "The total time for program with extraction and COA_CPU is " << total_Time_With_Extraction_CPU << endl;
	cout << "The total time for program without extraction and COA_CPU is " << total_Time_Without_Extraction_CPU << endl;
	cout << "The total time for program with extraction and COA_GPU is " << total_Time_With_Extraction_GPU<< endl;
	cout << "The total time for program without extraction and COA_GPU is " << total_Time_Without_Extraction_GPU<< endl;


	cout << "Press any key to exit" << endl;

	delete[] angleFileName;
	delete[] powerFileName;
	delete[] PMU_Angle;
	delete[] PMU_Power;

	//cout << "The total timetaken for the program to execute is " << total_GPU_Time << endl;

	cin.get();
	return 0;
}