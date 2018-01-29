#include <iostream>
#include <stdio.h>
#include <conio.h>
#include <string>
#include <fstream>
#include "pm.h"
#include <cmath>
#include <math.h>


using namespace std;

#define NUM_BUSES 11
#define SAMPLES 120
#define TOTAL_PMU 363
#define TIMESTEP 1/SAMPLES
#define VECTOR_SIZE_IN_BYTES_PMU (TOTAL_PMU*sizeof(float))
#define VECTOR_SIZE_IN_BYTES (TOTAL_PMU*SAMPLES*sizeof(float))


//Function to Generate file names
void GenerateFileNames(string *angle, string *power);

//Function to check the generated file names
void PrintFileNames(string *angle, string *power);

//Function to extract data from files
void ExtractData(string *angle, string *power, float *angleData, float *powerData);

//Function to extract line for a text file
float DataExtraction(char *data);

//Function to extract Center Of Angles
float CenterOfAnglesGPU(Performance& pm,float *angle, float *power, float *CenterOfAngle, float *product_PMU);

//Function to extract venter of angle
void CenterOfAngle(float *angle, float *power,float *COA);

//Function to predict stability
float StabilityGPU(Performance& pm,float *angle, float *COA);

