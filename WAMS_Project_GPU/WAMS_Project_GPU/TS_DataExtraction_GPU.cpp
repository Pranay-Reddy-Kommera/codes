#include "TS_Data_GPU.h"

void GenerateFileNames(string *angle, string *power)
{
	for (int i = 0; i < NUM_BUSES; i++)
	{
		angle[i] = string("Angle_Bus_") + to_string(i) + string(".txt");
		power[i] = string("RealPower_Bus_") + to_string(i) + string(".txt");
	}
}

void PrintFileNames(string *angle, string *power)
{
	for (int i = 0; i < NUM_BUSES; i++)
	{
		cout << angle[i];
		cout << "\n";
		cout << power[i];
		cout << "\n";
	}
}

void ExtractData(string *angle, string *power, float *angleData, float *powerData)
{
	char line[15];
	fstream input; //fstream object to read a file.

	//Angles Extraction
	for (int j = 0; j < NUM_BUSES; j++)
	{
		input.open(angle[j]);

		if (!input.is_open())
		{
			cout << angle[j] << " could not be opened" << endl;
			cout << "Press any key to exit";
			_getche();
		}

		for (int i = 0; i < SAMPLES; i++)
		{
			input.getline(line, 15);
			angleData[(i*TOTAL_PMU) + (j)] = DataExtraction(line);
		}

		input.close();

	}

	//Power Extraction
	for (int j = 0; j < NUM_BUSES; j++)
	{
		input.open(power[j]);

		if (!input.is_open())
		{
			cout << power[j] << " could not be opened" << endl;
			cout << "Press any key to exit";
			_getche();
		}

		for (int i = 0; i < SAMPLES; i++)
		{
			input.getline(line, 15);
			powerData[(i*TOTAL_PMU) + (j)] = DataExtraction(line);
		}

		input.close();

	}

	//Replicate data into all PMUs
	for (int n = 0; n < SAMPLES; n++)
	{
		for (int k = 0; k < ((TOTAL_PMU / NUM_BUSES) - 1); k++)
		{
			for (int l = 0; l < NUM_BUSES; l++)
			{
				angleData[(NUM_BUSES*(k + 1)) + l + (n*TOTAL_PMU)] = angleData[l + (n*TOTAL_PMU)];
				powerData[(NUM_BUSES*(k + 1)) + l + (n*TOTAL_PMU)] = powerData[l + (n*TOTAL_PMU)];

			}
		}
	}


}

float DataExtraction(char *data)
{
	char temp[20] = "";
	int k = 0;
	float derived_Data;

	//Read data
	int column_count = 0;
	while (column_count < 15)
	{
		if (data[column_count] != ' ')
		{
			temp[k] = data[column_count];
			k++;
		}
		column_count++;
	}
	temp[k] = '\0';
	derived_Data = strtof(temp, NULL);
	return derived_Data;
}

void CenterOfAngle(float *angle, float *power, float *COA)
{
	float sumPower = 0.0, productAnglePower = 0.0;

	for (int i = 0; i < TOTAL_PMU; i++)
	{
		sumPower += power[i];
		productAnglePower += (angle[i] * power[i]);
	}
	COA[0] = productAnglePower / sumPower;
}