#include "BA_data.h"
#include "BA_includeFiles.h"

__global__ void invVcuda(dtype* V_d, int num_pt)
{
	unsigned int k = threadIdx.x + (blockDim.x * threadIdx.y) + (blockDim.x * blockDim.y * blockIdx.x);
	
	if (k <= num_pt){
		dtype a, b, c, d, e, f, g, h, i, det; int Id = k * 9;
		a = V_d[Id]; b = V_d[Id + 1]; c = V_d[Id + 2];
		d = V_d[Id + 3]; e = V_d[Id + 4]; f = V_d[Id + 5];
		g = V_d[Id + 6]; h = V_d[Id + 7]; i = V_d[Id + 8];
		
		det = (a*((e*i) - (h*f))) - (b*((d*i) - (g*f))) + (c*((d*h) - (e*g)));

		V_d[Id] = ((e*i) - (h*f)) / det; V_d[Id + 1] = ((h*c) - (b*i)) / det; V_d[Id + 2] = ((b*f) - (e*c)) / det;
		V_d[Id + 3] = ((f*g) - (d*i)) / det; V_d[Id + 4] = ((a*i) - (g*c)) / det; V_d[Id + 5] = ((d*c) - (a*f)) / det;
		V_d[Id + 6] = ((d*h) - (e*g)) / det; V_d[Id + 7] = ((b*g) - (a*h)) / det; V_d[Id + 8] = ((a*e) - (b*d)) / det;
	}
}

//Kernel
void inverseVCalculation(dtype* V_d)
{
	//Call GPU kernel to perform V inverse
	int dimx = 32; int dimy = 32;
	cout << "Number of blocks is " << (num_pt / (dimx*dimy)) + 1 << endl;
	dim3 block(dimx, dimy); dim3 grid(((num_pt / (block.x*block.y)) + 1));

	invVcuda <<<grid, block >>>(V_d, num_pt);
	cudaDeviceSynchronize();
}

__global__ void tmpCuda(int* rowPtrU_d, int* colIndU_d, int* rowPtrWVWT_d, int* colIndWVWT_d, int* rowPtrS_d, int num_cam, int* tmp)
{
	unsigned int k = threadIdx.x + (blockDim.x * threadIdx.y) + (blockDim.x * blockDim.y * blockIdx.x);
	unsigned int newvwt, u, up1, wvwt, wvwtp1;
	unsigned int count, nes;

	if (k < (num_cam * 15))
	{
		wvwt = rowPtrWVWT_d[k]; wvwtp1 = rowPtrWVWT_d[k + 1]; newvwt = wvwtp1 - wvwt;
		u = rowPtrU_d[k]; up1 = rowPtrU_d[k + 1]; //neu = up1 - u;
		nes = newvwt;

		for (int i = u; i < up1; i++){
			count = 0;
			for (int j = wvwt; j < wvwtp1; j++){
				if (colIndU_d[i] == colIndWVWT_d[j]){
					count = count + 1;
				}
			}
			if (count == 0){
				nes = nes + 1;
			}
		}
		tmp[k] = nes;
		//printf("values of %d is %d\n", k, nes);
	}
	
	
	__syncthreads();

	
}

__global__ void rowPtrCuda(int* tmp, int* rowPtrS_d, int num_cam)
{
	unsigned int k = threadIdx.x + (blockDim.x * threadIdx.y) + (blockDim.x * blockDim.y * blockIdx.x);
	unsigned int total;

	if (k < (num_cam * 15)){
		
		if (k == 0){
			rowPtrS_d[0] = 0;
		}

		total = 0;
		for (int i = k; i >= 0; i--){
			total = total + tmp[i];
		}
		//printf("values of %d is %d\n", k, total);
		rowPtrS_d[k + 1] = total;
		//printf("values of %d is %d\n", k, rowPtrS_d[k+1]);
	}
	__syncthreads();
}

__global__ void scalCuda(int* rowPtrU_d, int* colIndU_d, dtype* U_d, int* rowPtrWVWT_d, int* colIndWVWT_d, dtype* WVWT_d, int* rowPtrS_d, int* colIndS_d, dtype* S_d, int num_cam)
{
	unsigned int k = threadIdx.x + (blockDim.x * threadIdx.y) + (blockDim.x * blockDim.y * blockIdx.x);
	unsigned int u, up1, wvwt, wvwtp1, s;
	unsigned int newvwt, num = 0;

	if (k < (num_cam * 15)){
		u = rowPtrU_d[k]; up1 = rowPtrU_d[k + 1]; //neu = up1 - u;
		wvwt = rowPtrWVWT_d[k]; wvwtp1 = rowPtrWVWT_d[k + 1]; newvwt = wvwtp1 - wvwt;
		s = rowPtrS_d[k]; //sp1 = rowPtrS_d[k + 1]; //nes = sp1 - s;

		for (int i = u; i < up1; i++){
			for (int j = (wvwt+num); j < wvwtp1; j++){
				if (colIndU_d[i] < colIndWVWT_d[j]){
					colIndS_d[s] = colIndU_d[i];
					S_d[s] = U_d[i];
					s = s + 1;
					break;
				}
				if (colIndU_d[i] == colIndWVWT_d[j]){
					colIndS_d[s] = colIndU_d[i];
					S_d[s] = U_d[i] - WVWT_d[j];
					s = s + 1;
					num = num + 1;
					break;
				}
				if (colIndU_d[i] > colIndWVWT_d[j]){
					colIndS_d[s] = colIndWVWT_d[j];
					S_d[s] = -WVWT_d[j];
					s = s + 1;
					num = num + 1;
				}
			}
		}

		if (num < newvwt){
			for (int i = (wvwt + num); i < wvwtp1; i++){
				colIndS_d[s] = colIndWVWT_d[i];
				S_d[s] = -WVWT_d[i];
				s = s + 1;
			}
		}
		//printf("values of %d is %d\n", k, rowPtrS_d[k + 1]);
	}
	__syncthreads();
}

//__global__ void printCuda(int* row, int size)
__global__ void printCuda(dtype* S_d, int size)
{
	for (int i = 0; i < 100; i++){
		printf("values of %d is %f\n", i, S_d[i]);
	}
}

void sCalculation(dtype* U_d, int* rowPtrU_d, int* colIndU_d, dtype* WVWT_d, int* rowPtrWVWT_d, int* colIndWVWT_d)
{
	cudaError_t cudaStat1, cudaStat2, cudaStat3, cudaStat4;
	//First fill out the rowPtrS_d as its size is row(U) + 1
	int* rowPtrS_d, *colIndS_d, *tmp; dtype* S_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrS_d, ((num_cam * 15) + 1)*sizeof(rowPtrS_d[0]));
	cudaStat2 = cudaMalloc((void**)&tmp, ((num_cam * 15))*sizeof(tmp[0]));
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess)){
		cout << "Successfully initialized rowPtrS_d and tmp" << endl;
	}
	else{
		cout << "****** Initialization of tmp and rowPtrS_d failed" << endl;
	}

	//Call GPU Kernel to fill rowPtrS_d
	int dimx = 16; int dimy = 16;
	cout << "Number of blocks is " << ((num_cam * 15) / (dimx*dimy)) + 1 << endl;
	dim3 block(dimx, dimy); dim3 grid((((num_cam * 15) / (block.x*block.y)) + 1));

	//rowPtrs_d is calculated
	tmpCuda << <grid, block >> >(rowPtrU_d, colIndU_d, rowPtrWVWT_d, colIndWVWT_d, rowPtrS_d, num_cam, tmp);
	cudaDeviceSynchronize();
	rowPtrCuda << <grid, block >> >(tmp, rowPtrS_d, num_cam);
	cudaDeviceSynchronize();
	cudaFree(tmp);

	//Allocate Memory for colIndS_d and S_d
	int jSId, baseS;
	cudaMemcpy(&jSId, rowPtrS_d + (num_cam * 15), sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&baseS, rowPtrS_d, sizeof(int), cudaMemcpyDeviceToHost);
	jSId -= baseS;

	//cout << "Value of jSId is " << jSId << " baseS is " << baseS << endl;
	cout << "The values of nnzS is " << jSId << " and total are " << (num_cam * num_cam * 225) << endl;

	cudaStat3 = cudaMalloc((void**)&colIndS_d, jSId*sizeof(colIndS_d[0]));
	cudaStat4 = cudaMalloc((void**)&S_d, jSId*sizeof(S_d[0]));
	if ((cudaStat3 == cudaSuccess) && (cudaStat4 == cudaSuccess)){
		cout << "Successfully initialized S_d and colIndS_d" << endl;
	}
	else{
		cout << "****** Initialization of S_d and colIndS_d failed" << endl;
	}

	//Calculate U - WVWT
	//Include U in S
	scalCuda << <grid, block >> >(rowPtrU_d, colIndU_d, U_d, rowPtrWVWT_d, colIndWVWT_d, WVWT_d, rowPtrS_d, colIndS_d, S_d,num_cam);
	cudaDeviceSynchronize();
	//dim3 block1(1, 1, 1), grid1(1, 1, 1);
	//printCuda << <grid1, block1 >> >(rowPtrS_d,((num_cam*15)+1));
	//printCuda << <grid1, block1 >> >(colIndS_d, jSId);
	//printCuda << <grid1, block1 >> >(S_d, jSId);
	//cudaDeviceSynchronize();

	cudaFree(U_d);
	cudaFree(rowPtrU_d);
	cudaFree(colIndU_d);
	cudaFree(WVWT_d);
	cudaFree(rowPtrWVWT_d);
	cudaFree(colIndWVWT_d);

	//Allocate CPU memory for S
	dtype* S = new dtype[jSId];
	int* rowPtrS = new int[(num_cam * 15) + 1];
	int* colIndS = new int[jSId];

	cudaStat1 = cudaMemcpy(S, S_d, jSId*sizeof(dtype), cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(rowPtrS, rowPtrS_d, ((num_cam * 15) + 1)*sizeof(int), cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(colIndS, colIndS_d, jSId*sizeof(int), cudaMemcpyDeviceToHost);
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess)){
		cout << "Successfully copied S into host" << endl;
	}
	else{
		cout << "****** Copy of S into host failed" << endl;
	}

	cudaFree(S_d);
	cudaFree(rowPtrS_d);
	cudaFree(colIndS_d);

	//Print S and verify
	/*for (int i = 0; i < 100; i++){
		cout << rowPtrS[i] << "\t";
	}*/
	cout << endl;
	cout << "************************************************************************" << endl;
	////////////////////////////////////////////////////////////////////////////////////////////////
}