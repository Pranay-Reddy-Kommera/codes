#include <cuda.h>
#include "TS_Data_GPU.h"


bool HandleCUDAErrors(cudaError_t t)
{
	if (t != cudaSuccess)
	{
		puts(cudaGetErrorString(cudaGetLastError()));
		return false;
	}
	return true;
}

bool GetCUDARunTimeError()
{
	cudaError_t t = cudaGetLastError();
	if (t != cudaSuccess)
	{
		puts(cudaGetErrorString(t));
		return false;
	}
	return true;
}

__global__ void ProductPMU(float* angle, float* power, float* prod)
{
	unsigned int x = threadIdx.x + (threadIdx.y*blockDim.x);
	unsigned int y = (blockDim.x*blockDim.y*blockIdx.x) + x;

	prod[y] = angle[y] * power[y];

}

float CenterOfAnglesGPU(Performance& pm, float *angle, float *power, float *CenterOfAngle, float *product_PMU)
{
	
	float *angleGPU, *powerGPU, *productGPU, add=0.0, powerAdd=0.0,elapsedTime=0.0;

	//Allocate memory on GPU
	HandleCUDAErrors(cudaMalloc((void**)&angleGPU, VECTOR_SIZE_IN_BYTES_PMU));
	HandleCUDAErrors(cudaMalloc((void**)&powerGPU, VECTOR_SIZE_IN_BYTES_PMU));
	HandleCUDAErrors(cudaMalloc((void**)&productGPU, VECTOR_SIZE_IN_BYTES_PMU));

	pm.OnSetStartTime();
	//Copy data into device memory
	HandleCUDAErrors(cudaMemcpy(angleGPU, angle, VECTOR_SIZE_IN_BYTES_PMU, cudaMemcpyHostToDevice));
	HandleCUDAErrors(cudaMemcpy(powerGPU, power, VECTOR_SIZE_IN_BYTES_PMU, cudaMemcpyHostToDevice));
	

	//Grid and Block Setup
	int blockx = 16;
	int blocky = 16;
	dim3 block(blockx, blocky);
	dim3 grid((TOTAL_PMU+255)/256, 1);

	pm.OnSetStopTime();
	float time_Copy = pm.OnGetElapsedTime();
	//cout << "The copy time is " << time_Copy << endl;

	cudaEvent_t kernel_start; //Set the event corresponding to kernel lanuch
	cudaEvent_t kernel_stop; //Set the event corresponding to kernel completion

	HandleCUDAErrors(cudaEventCreate(&kernel_start));//start event object is created here
	HandleCUDAErrors(cudaEventCreate(&kernel_stop));//stop event object is created here

	HandleCUDAErrors(cudaEventRecord(kernel_start, 0));//recording the kernel start, 0 is a string//its a non blocking call

	//Kernel Calling
	ProductPMU << <grid, block >> >(angleGPU, powerGPU, productGPU);
	//cudaThreadSynchronize();
	GetCUDARunTimeError();

	HandleCUDAErrors(cudaEventRecord(kernel_stop, 0));//record the kernel stop
	HandleCUDAErrors(cudaEventSynchronize(kernel_stop));//This is a blocking call same as cudaThreadSynchronize

	//Used to store elapsed time
	HandleCUDAErrors(cudaEventElapsedTime(&elapsedTime, kernel_start, kernel_stop));//address is where we store the value of elapsed time

	//cout << "The kernel time is" << elapsedTime << endl;
	pm.OnSetStartTime();
	HandleCUDAErrors(cudaMemcpy(product_PMU, productGPU, VECTOR_SIZE_IN_BYTES_PMU, cudaMemcpyDeviceToHost));

	for (int i = 0; i < TOTAL_PMU; i++)
	{
		add += product_PMU[i];
		powerAdd += power[i];
	}

	CenterOfAngle[0] = add / powerAdd;
	pm.OnSetStopTime();
	float time_CopyBack = pm.OnGetElapsedTime();
	//cout << "The copy back time is " << time_CopyBack << endl;

	float time_Total = (time_Copy * 1000) + elapsedTime + (time_CopyBack * 1000);
	//cout << "The total time is " << time_Total<< endl;
	//cout << "The center of angle is " << CenterOfAngle[0] << endl;

	HandleCUDAErrors(cudaFree(angleGPU));
	HandleCUDAErrors(cudaFree(powerGPU));
	HandleCUDAErrors(cudaFree(productGPU));

	//Remove the kernel events
	HandleCUDAErrors(cudaEventDestroy(kernel_start));
	HandleCUDAErrors(cudaEventDestroy(kernel_stop));
	HandleCUDAErrors(cudaDeviceReset());

	return time_Total;
}


__global__ void Stability(float* d_Angle, float positiveThreshold, float negativeThreshold)
{
	unsigned int x = threadIdx.x + (threadIdx.y*blockDim.x);
	unsigned int y = (blockDim.x*blockDim.y*blockIdx.x) + x;
	//float neg = negativeThreshold;
	//float pos = positiveThreshold;
	
		int temp = 0;
		float area = 0.0;

		for (int i = 0; i < (SAMPLES - 1); i++)
		{
			float anglePrevious = d_Angle[(i*TOTAL_PMU) + y];
			float angle = d_Angle[((i+1)*TOTAL_PMU)+y];
			float diffAngles = angle - anglePrevious;


			if ((angle <= positiveThreshold) && (angle >= negativeThreshold))
			{
				temp += 0;
			}
			else
			{
				if (angle > positiveThreshold)
				{
					if (anglePrevious <= positiveThreshold)
					{
						anglePrevious = positiveThreshold;
					}

					//area += (fabsf(anglePrevious)*TIMESTEP) + (0.5*TIMESTEP*fabsf(diffAngles));

					if (diffAngles < 0)
					{
					diffAngles = -diffAngles;
					}
					if (anglePrevious < 0)
					{
					anglePrevious = -anglePrevious;
					}
					area += (anglePrevious*TIMESTEP) + (0.5*TIMESTEP*(diffAngles));
					if (area >= 5)
					{
						temp += 1;
					}
					else
					{
						temp += 0;
						area = 0.0;
					}
				}
				else
				{
					if (anglePrevious >= negativeThreshold)
					{
						anglePrevious = negativeThreshold;
					}

					//area += (fabsf(anglePrevious)*TIMESTEP) + (0.5*TIMESTEP*fabsf(diffAngles));

					if (diffAngles < 0)
					{
					diffAngles = -diffAngles;
					}
					if (anglePrevious < 0)
					{
					anglePrevious = -anglePrevious;
					}
					area += (anglePrevious*TIMESTEP) + (0.5*TIMESTEP*(diffAngles));
					if (area >= 5)
					{
						temp += 1;
					}
					else
					{
						temp += 0;
						area = 0.0;
					}
				}
			}
		}
		if ((temp > 0)&&(y<TOTAL_PMU))
		{
			//cout << "The system is approaching unstable consition at generator bus number " << j << endl;
			printf("The system is approaching unstable condition at generator bus number %d\n", y);
		}
	
}

float StabilityGPU(Performance& pm, float *angle, float *COA)
{
	float *d_Angle, elapsedTime=0.0;

	//Allocate device memory
	HandleCUDAErrors(cudaMalloc((void**)&d_Angle, VECTOR_SIZE_IN_BYTES));

	//Copy data to GPU
	pm.OnSetStartTime();
	HandleCUDAErrors(cudaMemcpy(d_Angle, angle, VECTOR_SIZE_IN_BYTES, cudaMemcpyHostToDevice));
	

	//Grid and Block Setup
	int blockx = 16;
	int blocky = 16;
	dim3 block(blockx, blocky);
	dim3 grid((TOTAL_PMU + 255) / 256, 1);

	float positiveThreshold = COA[0] + 60;
	float negativeThreshold = COA[0] - 70;

	pm.OnSetStopTime();
	float time_Copy_Stability = pm.OnGetElapsedTime();
	//cout << "The copy time is " << time_Copy_Stability << endl;

	//cout << "The +ve angle is " << positiveThreshold << endl;
	//cout << "The -ve threshold is " << negativeThreshold << endl;
	
	cudaEvent_t kernel_start; //Set the event corresponding to kernel lanuch
	cudaEvent_t kernel_stop; //Set the event corresponding to kernel completion

	HandleCUDAErrors(cudaEventCreate(&kernel_start));//start event object is created here
	HandleCUDAErrors(cudaEventCreate(&kernel_stop));//stop event object is created here

	HandleCUDAErrors(cudaEventRecord(kernel_start, 0));//recording the kernel start, 0 is a string//its a non blocking call

	Stability << <grid, block >> >(d_Angle,positiveThreshold, negativeThreshold);
	//cudaThreadSynchronize();
	GetCUDARunTimeError();

	HandleCUDAErrors(cudaEventRecord(kernel_stop, 0));//record the kernel stop
	HandleCUDAErrors(cudaEventSynchronize(kernel_stop));//This is a blocking call same as cudaThreadSynchronize

	//Used to store elapsed time
	HandleCUDAErrors(cudaEventElapsedTime(&elapsedTime, kernel_start, kernel_stop));//address is where we store the value of elapsed time
	//cout << "The elapsed time is " << elapsedTime << endl;
	float total_Time = (time_Copy_Stability * 1000) + elapsedTime;
	HandleCUDAErrors(cudaFree(d_Angle));

	//Remove the kernel events
	HandleCUDAErrors(cudaEventDestroy(kernel_start));
	HandleCUDAErrors(cudaEventDestroy(kernel_stop));
	HandleCUDAErrors(cudaDeviceReset());

	return total_Time;
}
