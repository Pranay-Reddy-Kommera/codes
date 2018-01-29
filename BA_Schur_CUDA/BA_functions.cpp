///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_functions.cpp
///   Description    :         Functions required in Bundle Adjustment
//////////////////////////////////////////////////////////////////////////////////////

#include "BA_includeFiles.h"
#include "BA_data.h"

//Calculate Measured projections 
//Calculating f(P)
void MeasuredProjections()
{
	int camId, ptId, projId;
	dtype rxpt1, rxpt2, rxpt3, pDiv[2], ppNorm, frp;
	//Calculating the measured projections
	for (int i = 0; i < num_proj; i++){

		//Calculate RX + T
		camId = (camidx[i]) * 15; ptId = (num_cam * 15) + ((ptidx[i]) * 3);
		rxpt1 = (P[ptId] * P[camId]) + (P[ptId + 1] * P[camId + 1]) + (P[ptId + 2] * P[camId + 2]) + P[camId + 9];
		rxpt2 = (P[ptId] * P[camId + 3]) + (P[ptId + 1] * P[camId + 4]) + (P[ptId + 2] * P[camId + 5]) + P[camId + 10];
		rxpt3 = (P[ptId] * P[camId + 6]) + (P[ptId + 1] * P[camId + 7]) + (P[ptId + 2] * P[camId + 8]) + P[camId + 11];

		//Calculate Perspective Division; p = -P/P.z
		pDiv[0] = -rxpt1 / rxpt3; pDiv[1] = -rxpt2 / rxpt3;

		//Calculating pixel coordinates; p'=f * r(p) * p
		ppNorm = vecNorm(pDiv, 2);
		frp = P[camId + 12] * (1 + ((P[camId + 13])*(ppNorm*ppNorm)) + ((P[camId + 14])*(ppNorm*ppNorm*ppNorm*ppNorm)));

		//Store in calcProjections
		projId = 2 * i;
		calcProjections[projId] = frp*pDiv[0]; calcProjections[projId + 1] = frp*pDiv[1];
	}
}

void MeasuredProjectionsnew()
{
	int camId, ptId, projId;
	dtype rxpt1, rxpt2, rxpt3, pDiv[2], ppNorm, frp;
	//Calculating the measured projections
	for (int i = 0; i < num_proj; i++){

		//Calculate RX + T
		camId = (camidx[i]) * 15; ptId = (num_cam * 15) + ((ptidx[i]) * 3);
		rxpt1 = (Pnew[ptId] * Pnew[camId]) + (Pnew[ptId + 1] * Pnew[camId + 1]) + (Pnew[ptId + 2] * Pnew[camId + 2]) + Pnew[camId + 9];
		rxpt2 = (Pnew[ptId] * Pnew[camId + 3]) + (Pnew[ptId + 1] * Pnew[camId + 4]) + (Pnew[ptId + 2] * Pnew[camId + 5]) + Pnew[camId + 10];
		rxpt3 = (Pnew[ptId] * Pnew[camId + 6]) + (Pnew[ptId + 1] * Pnew[camId + 7]) + (Pnew[ptId + 2] * Pnew[camId + 8]) + Pnew[camId + 11];

		//Calculate Perspective Division; p = -P/P.z
		pDiv[0] = -rxpt1 / rxpt3; pDiv[1] = -rxpt2 / rxpt3;

		//Calculating pixel coordinates; p'=f * r(p) * p
		ppNorm = vecNorm(pDiv, 2);
		frp = Pnew[camId + 12] * (1 + ((Pnew[camId + 13])*(ppNorm*ppNorm)) + ((Pnew[camId + 14])*(ppNorm*ppNorm*ppNorm*ppNorm)));

		//Store in calcProjections
		projId = 2 * i;
		calcProjectionsnew[projId] = frp*pDiv[0]; calcProjectionsnew[projId + 1] = frp*pDiv[1];
	}
}

//Calculate RMS erro values
void errorCalculation()
{
	for (int i = 0; i < (2 * num_proj); i++){
		Ep[i] = projections[i] - calcProjections[i];
	}
}

void RMSErrorCalculation()
{
	dtype EpNorm, tmp, rms;

	EpNorm = vecNorm(Ep, 2 * num_proj);
	tmp = (EpNorm*EpNorm) / (2 * num_proj);
	rms = sqrt(tmp);
	errorVector.push_back(rms);

	//cout << "RMS Error Norm is " << rms << endl;

}

//Calculate U, V, W, Ea, Eb
void UVWgCalculation()
{
	cudaError_t cudaStat1, cudaStat2, cudaStat3, cudaStat4, cudaStat5, cudaStat6;
	cusparseStatus_t status, status1, status2, status3, status4, status5, status6;
	cusparseHandle_t handle = 0;
	cusparseMatDescr_t descr = 0;

	long long int jcId = num_proj * 2 * 15; long long int jpId = num_proj * 2 * 3; long long int num_row = (2 * num_cam*num_pt); long long int num_colJc = 15 * num_cam; long long int num_colJp = 3 * num_pt;
	int dsize = sizeof(dtype);
	int isize = sizeof(int);
	//Allocate memory for Jacobian
	//Jacobian are first stored in COO format
	dtype* Jc = new dtype[jcId];
	dtype* Jp = new dtype[jpId];

	//COO Row and Column Indices
	int* rowJc = new int[jcId];
	int* colJc = new int[jcId];
	int* rowJp = new int[jpId];
	int* colJp = new int[jpId];

	//Cuda Initialization
	int* rowJc_d, *rowJp_d, *colIndJc_d, *colIndJp_d; dtype* Jc_d, *Jp_d;
	int* rowPtrJc_d, *rowPtrJp_d;

	//Calculate Jc and Jp
	NormalEquations(Jc, Jp, rowJc, colJc, rowJp, colJp);

	//Checking 
	/*for(int i = 0; i < 200; i++){
	cout << rowJc[i] << "\t";
	cout << endl;
	}*/
	cout << "************************************************************************" << endl;
	cout << "************************************************************************" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// COO cuda Allocation
	cudaStat1 = cudaMalloc((void**)&rowJc_d, jcId*sizeof(rowJc_d[0]));
	cudaStat2 = cudaMalloc((void**)&rowJp_d, jpId*sizeof(rowJp_d[0]));
	cudaStat3 = cudaMalloc((void**)&colIndJc_d, jcId*sizeof(colIndJc_d[0]));
	cudaStat4 = cudaMalloc((void**)&colIndJp_d, jpId*sizeof(colIndJp_d[0]));
	cudaStat5 = cudaMalloc((void**)&Jc_d, jcId*sizeof(Jc_d[0]));
	cudaStat6 = cudaMalloc((void**)&Jp_d, jpId*sizeof(Jp_d[0]));

	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat4 == cudaSuccess) && (cudaStat5 == cudaSuccess) && (cudaStat6 == cudaSuccess)){
	//if ((cudaStat1 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat5 == cudaSuccess)){
		cout << "1.Successfully Allocated COO in GPU memory" << endl;
	}
	else{
		cout << "*** No (1) *** Did not allocate all COO variables in GPU Memory" << endl;
	}

	int total = ((2 * jcId*isize) + (2 * jpId*isize) + ((jcId + jpId)*dsize)) / (1024 * 1024);
	cout << "001.Total COO memory Allocated is " << total << " MB" << endl;

	//CSR cuda Allocation
	cudaStat1 = cudaMalloc((void**)&rowPtrJc_d, (num_row + 1)*sizeof(rowPtrJc_d[0]));
	cudaStat2 = cudaMalloc((void**)&rowPtrJp_d, (num_row + 1)*sizeof(rowPtrJp_d[0]));

	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess)){
	//if ((cudaStat1 == cudaSuccess)){
		cout << "2.Successfully Allocated CSR rowPtr in GPU memory" << endl;
	}
	else{
		cout << "*** No (2) *** Did not allocate all CSR variables in GPU Memory" << endl;
	}

	total = total + ((2 * (num_row + 1)*isize) / (1024 * 1024));
	cout << "002.Total COO + rowPtrJc and Jp memory Allocated is " << total << " MB" << endl;

	//Copy CPU data into GPU
	cudaStat1 = cudaMemcpy(rowJc_d, rowJc, (jcId*sizeof(rowJc_d[0])), cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(rowJp_d, rowJp, (jpId*sizeof(rowJp_d[0])), cudaMemcpyHostToDevice);
	cudaStat3 = cudaMemcpy(colIndJc_d, colJc, (jcId*sizeof(colIndJc_d[0])), cudaMemcpyHostToDevice);
	cudaStat4 = cudaMemcpy(colIndJp_d, colJp, (jpId*sizeof(colIndJp_d[0])), cudaMemcpyHostToDevice);
	cudaStat5 = cudaMemcpy(Jc_d, Jc, (jcId*sizeof(Jc_d[0])), cudaMemcpyHostToDevice);
	cudaStat6 = cudaMemcpy(Jp_d, Jp, (jpId*sizeof(Jp_d[0])), cudaMemcpyHostToDevice);

	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat4 == cudaSuccess) && (cudaStat5 == cudaSuccess) && (cudaStat6 == cudaSuccess)){
	//if ((cudaStat1 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat5 == cudaSuccess)){
		cout << "3.Successfully copied COO into GPU memory" << endl;
	}
	else{
		cout << "*** No (3) *** Did not copy COO variables into GPU Memory" << endl;
	}

	//Initialize cusparse library
	status = cusparseCreate(&handle);
	if (status == CUSPARSE_STATUS_SUCCESS){
		cout << "4.Initialized cusparse successfully" << endl;
	}
	else{
		cout << "*** No (4) *** cusparse initialization failed" << endl;
	}

	//Create and setup matrix descriptor
	status1 = cusparseCreateMatDescr(&descr);
	if (status1 == CUSPARSE_STATUS_SUCCESS){
		cout << "5.Initialized matrix descriptor successfully" << endl;
	}
	else{
		cout << "*** No (5) *** Matrix descriptor initialization failed" << endl;
	}

	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

	//Convert COO to CSR
	status2 = cusparseXcoo2csr(handle, rowJc_d, jcId, num_row, rowPtrJc_d, CUSPARSE_INDEX_BASE_ZERO);
	status3 = cusparseXcoo2csr(handle, rowJp_d, jpId, num_row, rowPtrJp_d, CUSPARSE_INDEX_BASE_ZERO);
	cudaDeviceSynchronize();

	if ((status2 == CUSPARSE_STATUS_SUCCESS) && (status3 == CUSPARSE_STATUS_SUCCESS)){
	//if ((status1 == CUSPARSE_STATUS_SUCCESS)){
		cout << "6.COO to CSR conversion done successfully" << endl;
	}
	else{
		cout << "*** No (6) *** COO to CSR conversion failed" << endl;
	}
	cout << "************************************************************************" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Free COO elements which are not needed
	//Both on CPU and GPU
	cudaFree(rowJc_d);
	cudaFree(rowJp_d);
	delete[] rowJc;
	delete[] rowJp;
	delete[] colJc;
	delete[] colJp;
	delete[] Jc;
	delete[] Jp;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Initialize U and calculate
	//Initialize and Allocate rowPtrJT_c
	int* rowPtrU_d, *colIndU_d; dtype* U_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrU_d, (num_colJc + 1)*sizeof(rowPtrU_d[0]));

	if ((cudaStat1 == cudaSuccess)){
		cout << "7.Successfully Allocated CSR rowPtr for U in GPU memory" << endl;
	}
	else{
		cout << "*** No (7) *** CSR rowPtr for U allocation failed" << endl;
	}

	//Create Jc^T in CSR format. Jc^T in CSR == Jc in CSC
	//Initialize and Allocate csrJ^T
	dtype* JcT_d; int* rowPtrJcT_d, *colIndJcT_d;
	cudaStat2 = cudaMalloc((void**)&rowPtrJcT_d, (num_colJc + 1)*sizeof(rowPtrJcT_d[0]));
	cudaStat3 = cudaMalloc((void**)&JcT_d, jcId*sizeof(JcT_d[0]));
	cudaStat4 = cudaMalloc((void**)&colIndJcT_d, jcId*sizeof(colIndJcT_d[0]));
	if ((cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat4 == cudaSuccess)){
		cout << "8.Successfully initialized Jc^T in GPU memory" << endl;
	}
	else{
		cout << "*** No (8) *** Jc^T initialization in GPU Memory failed" << endl;
	}
	cudaDeviceSynchronize();
	status4 = cusparseDcsr2csc(handle, num_row, num_colJc, jcId, Jc_d, rowPtrJc_d, colIndJc_d, JcT_d, colIndJcT_d, rowPtrJcT_d, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
	if ((status4 == CUSPARSE_STATUS_SUCCESS)){
		cout << "9.CSR to CSC conversion done successfully for Jc" << endl;
	}
	else{
		cout << "*** No (9) *** CSR to CSC conversion failed" << endl;
	}

	cudaDeviceSynchronize();
	int jUId, baseU;
	int* nnzTotalDevHostPtrU = &jUId;
	status5 = cusparseXcsrgemmNnz(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJc, num_row, descr, jcId, rowPtrJcT_d, colIndJcT_d, descr, jcId, rowPtrJc_d, colIndJc_d, descr, rowPtrU_d, nnzTotalDevHostPtrU);
	if ((status5 == CUSPARSE_STATUS_SUCCESS)){
		cout << "10.NNz in U is done successfully" << endl;
	}
	else{
		cout << "*** No (10) *** NNz in U failed" << endl;
	}

	cudaDeviceSynchronize();
	if (NULL != nnzTotalDevHostPtrU){
		jUId = *nnzTotalDevHostPtrU;
	}
	else{
		cudaMemcpy(&jUId, rowPtrU_d + num_colJc, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseU, rowPtrU_d, sizeof(int), cudaMemcpyDeviceToHost);
		jUId -= baseU;
	}

	cout << "11.The values of nnzU is " << jUId << " and total are " << (num_cam * 225) << endl;
	if (jUId == (num_cam * 225)){
		cout << "12.All the U matrices are filled" << endl;
	}
	else{
		cout << "12.*** All U are not filled" << endl;
	}

	cudaStat5 = cudaMalloc((void**)&colIndU_d, jUId*sizeof(colIndU_d[0]));
	cudaStat6 = cudaMalloc((void**)&U_d, jUId*sizeof(U_d[0]));
	if ((cudaStat5 == cudaSuccess) && (cudaStat6 == cudaSuccess)){
		cout << "13.Successfully initialized U and colIndU" << endl;
	}
	else{
		cout << "*** No (13) *** Initialization of U and colIndU failed" << endl;
	}

	cudaDeviceSynchronize();
	//Finally calculating U
	status6 = cusparseDcsrgemm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJc, num_row, descr, jcId, JcT_d, rowPtrJcT_d, colIndJcT_d, descr, jcId, Jc_d, rowPtrJc_d, colIndJc_d, descr, U_d, rowPtrU_d, colIndU_d);
	if ((status6 == CUSPARSE_STATUS_SUCCESS)){
		cout << "14.U is calculated successfully" << endl;
	}
	else{
		cout << "*** No (14) *** U calculation failed" << endl;
	}
	cudaDeviceSynchronize();
	cout << "************************************************************************" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Free Jc data
	cudaFree(Jc_d);
	cudaFree(rowPtrJc_d);
	cudaFree(colIndJc_d);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Initialize V and calculate
	//Initialize and Allocate rowPtrJT_p
	int* rowPtrV_d, *colIndV_d; dtype* V_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrV_d, (num_colJp + 1)*sizeof(rowPtrV_d[0]));

	if ((cudaStat1 == cudaSuccess)){
		cout << "15.Successfully Allocated CSR rowPtr for V in GPU memory" << endl;
	}
	else{
		cout << "*** No (15) *** CSR rowPtr for V in GPU allocation failed" << endl;
	}

	//Create Jp^T in CSR format. Jp^T in CSR == Jp in CSC
	//Initialize and Allocate csrJ^T
	dtype* JpT_d; int* rowPtrJpT_d, *colIndJpT_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrJpT_d, (num_colJp + 1)*sizeof(rowPtrJpT_d[0]));
	cudaStat2 = cudaMalloc((void**)&JpT_d, jpId*sizeof(JpT_d[0]));
	cudaStat3 = cudaMalloc((void**)&colIndJpT_d, jpId*sizeof(colIndJpT_d[0]));
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess)){
		cout << "16.Successfully initialized Jp^T in GPU memory" << endl;
	}
	else{
		cout << "*** No (16) *** Initialization of Jp^T in GPU Memory failed" << endl;
	}

	cudaDeviceSynchronize();
	status2 = cusparseDcsr2csc(handle, num_row, num_colJp, jpId, Jp_d, rowPtrJp_d, colIndJp_d, JpT_d, colIndJpT_d, rowPtrJpT_d, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "17.CSR to CSC conversion done successfully for Jp" << endl;
	}
	else{
		cout << "*** No (17) *** CSR to CSC conversion for Jp failed" << endl;
	}

	cudaDeviceSynchronize();
	int jVId, baseV;
	int* nnzTotalDevHostPtr = &jVId;
	status2 = cusparseXcsrgemmNnz(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJp, num_colJp, num_row, descr, jpId, rowPtrJpT_d, colIndJpT_d, descr, jpId, rowPtrJp_d, colIndJp_d, descr, rowPtrV_d, nnzTotalDevHostPtr);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "18.NNz in V is done successfully" << endl;
	}
	else{
		cout << "*** No (18) *** NNz in V failed" << endl;
	}

	cudaDeviceSynchronize();
	if (NULL != nnzTotalDevHostPtr){
		jVId = *nnzTotalDevHostPtr;
	}
	else{
		cudaMemcpy(&jVId, rowPtrV_d + num_colJp, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseV, rowPtrV_d, sizeof(int), cudaMemcpyDeviceToHost);
		jVId -= baseV;
	}
	
	cout << "19.The values of nnzU is " << jVId << " and total are " << (num_pt * 9) << endl;
	if (jVId == (num_pt * 9)){
		cout << "20.All V matrices are filled" << endl;
	}
	else{
		cout << "20.*** All V are not filled" << endl;
	}

	cudaStat1 = cudaMalloc((void**)&colIndV_d, jVId*sizeof(colIndV_d[0]));
	cudaStat2 = cudaMalloc((void**)&V_d, jVId*sizeof(V_d[0]));
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess)){
		cout << "21.Successfully initialized V and colIndV" << endl;
	}
	else{
		cout << "*** No (21) *** Initialization of V and colIndV failed" << endl;
	}

	cudaDeviceSynchronize();
	//Finally calculating U
	status2 = cusparseDcsrgemm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJp, num_colJp, num_row, descr, jpId, JpT_d, rowPtrJpT_d, colIndJpT_d, descr, jpId, Jp_d, rowPtrJp_d, colIndJp_d, descr, V_d, rowPtrV_d, colIndV_d);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "22.V is calculated successfully" << endl;
	}
	else{
		cout << "*** No (22) *** V calculation failed" << endl;
	}

	cudaDeviceSynchronize();
	cout << "************************************************************************" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Free Jc data
	cudaFree(JpT_d);
	cudaFree(rowPtrJpT_d);
	cudaFree(colIndJpT_d);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Initialize W and calculate
	//Initialize and Allocate rowPtrJT_p
	int* rowPtrW_d, *colIndW_d; dtype* W_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrW_d, (num_colJc + 1)*sizeof(rowPtrW_d[0]));

	if ((cudaStat1 == cudaSuccess)){
		cout << "23.Successfully Allocated CSR rowPtr for W in GPU memory" << endl;
	}
	else{
		cout << "*** No (23) *** Allocation of CS rowPtr for W in GPU Memory failed" << endl;
	}

	cudaDeviceSynchronize();
	int jWId, baseW;
	int* nnzTotalDevHostPtrW = &jWId;
	status2 = cusparseXcsrgemmNnz(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJp, num_row, descr, jcId, rowPtrJcT_d, colIndJcT_d, descr, jpId, rowPtrJp_d, colIndJp_d, descr, rowPtrW_d, nnzTotalDevHostPtrW);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "24.NNz in W is done successfully" << endl;
	}
	else{
		cout << "*** No (24) *** NNz in W failed" << endl;
	}

	cudaDeviceSynchronize();
	if (NULL != nnzTotalDevHostPtrW){
		jWId = *nnzTotalDevHostPtrW;
	}
	else{
		cudaMemcpy(&jWId, rowPtrW_d + num_colJc, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseW, rowPtrW_d, sizeof(int), cudaMemcpyDeviceToHost);
		jWId -= baseW;
	}
	
	cout << "25.The values of nnzW is " << jWId << " and total are " << ((15 * num_cam) * (3 * num_pt)) << endl;
	if (jWId == ((15 * num_cam) * (3 * num_pt))){
		cout << "26.All W elements are filled" << endl;
	}
	else{
		cout << "26.*** All W are not filled" << endl;
	}

	cudaStat1 = cudaMalloc((void**)&colIndW_d, jWId*sizeof(colIndW_d[0]));
	cudaStat2 = cudaMalloc((void**)&W_d, jWId*sizeof(W_d[0]));
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess)){
		cout << "27.Successfully initialized W and colIndW" << endl;
	}
	else{
		cout << "*** No (27) *** Initialization of W and colIndW failed" << endl;
	}

	cudaDeviceSynchronize();
	//Finally calculating U
	status2 = cusparseDcsrgemm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJp, num_row, descr, jcId, JcT_d, rowPtrJcT_d, colIndJcT_d, descr, jpId, Jp_d, rowPtrJp_d, colIndJp_d, descr, W_d, rowPtrW_d, colIndW_d);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "28.W is calculated successfully" << endl;
	}
	else{
		cout << "*** No (28) *** W calculation failed" << endl;
	}

	cudaDeviceSynchronize();
	cout << "************************************************************************" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Free Jc data
	cudaFree(Jp_d);
	cudaFree(rowPtrJp_d);
	cudaFree(colIndJp_d);
	cudaFree(JcT_d);
	cudaFree(rowPtrJcT_d);
	cudaFree(colIndJcT_d);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Initialize U, V and W on CPU
	int* rowPtrU = new int[num_colJc + 1];
	int* rowPtrV = new int[num_colJp + 1];
	//int* colIndU = new int[jUId];
	//dtype* W = new dtype[jWId];

	//Copy U details to host
	cudaStat1 = cudaMemcpy(U, U_d, jUId*sizeof(U_d[0]), cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(V, V_d, jVId*sizeof(V_d[0]), cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(rowPtrU, rowPtrU_d, (num_colJc + 1)*sizeof(int), cudaMemcpyDeviceToHost);
	cudaStat4 = cudaMemcpy(rowPtrV, rowPtrV_d, (num_colJp + 1)*sizeof(int), cudaMemcpyDeviceToHost);
	//cudaStat3 = cudaMemcpy(colIndU, colIndU_d, jUId*sizeof(int), cudaMemcpyDeviceToHost);
	//cudaStat2 = cudaMemcpy(V, V_d, jVId*sizeof(V_d[0]), cudaMemcpyDeviceToHost);
	//cudaStat3 = cudaMemcpy(W, W_d, jWId*sizeof(W_d[0]), cudaMemcpyDeviceToHost);

	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat4 == cudaSuccess)){
		cout << "29.U, V and RowPtr of U and V are copied onto host" << endl;
	}
	else{
		cout << "*** No (29) *** U, V and rowPtr of U and V are not copied to host" << endl;
	}

	int count = 0;
	for (int i = 0; i < (num_colJc); i++){
		if ((rowPtrU[i + 1] - rowPtrU[i]) != 15){
			count = count + 1;
		}
	}
	if (count == 0){
		cout << "30.U has expected layout" << endl;
	}
	else{
		cout << "30.U has unexpected layout" << endl;
	}
	count = 0;
	for (int i = 0; i < (num_colJp); i++){
		if ((rowPtrV[i + 1] - rowPtrV[i]) != 3){
			count = count + 1;
		}
	}
	if (count == 0){
		cout << "31.V has expected layout" << endl;
	}
	else{
		cout << "31.V has unexpected layout" << endl;
	}
	//delete rowPtr
	delete[] rowPtrU;
	delete[] rowPtrV;

	cout << endl;
	cout << "************************************************************************" << endl;
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Total memory U, V and W occupies in GPU is 
	int valSize = (jUId*dsize) + (jVId*dsize) + (jWId*dsize);
	int colIndSize = (jUId*isize) + (jVId*isize) + (jWId*isize);
	int rowPtrSize = ((num_colJc + 1)*isize) + ((num_colJp + 1)*isize) + ((num_colJc + 1)*isize);

	cout << "32.Total memory utilized for U, V and W is " << ((valSize + colIndSize + rowPtrSize) / (1024 * 1024)) << " MB" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

auto start = chrono::high_resolution_clock::now();

	//Calculate mu
	mu = muCalculation();
	cout << "mu is " << mu << endl;

	//Calculate augmented U and V matrix
	augmentUV();

	//invVCalculation();

	//Copy U and V back to GPU
	cudaStat1 = cudaMemcpy(U_d, U, (jUId*sizeof(U_d[0])), cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(V_d, V, (jVId*sizeof(U_d[0])), cudaMemcpyHostToDevice);
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess)){
		cout << "33.Successfully copied U and V onto GPU Memory" << endl;
	}
	else{
		cout << "*** No (33) *** Copy of U and V onto GPU Memory failed" << endl;
	}

	inverseVCalculation(V_d);

	cout << "mu, augmentation and inverse of V are calculated" << endl;
	/*
	//Copy Vinv to CPU
	cudaStat1 = cudaMemcpy(V, V_d, jVId*sizeof(V_d[0]), cudaMemcpyDeviceToHost);

	for (int i = 0; i < 100; i++){
		cout << V[i] << "\t";
	}
	cout << endl;
	*/
	cout << "************************************************************************" << endl;

	///////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate S = U - (WVW')
	//First calculate W * V
	//Initialize and Allocate rowPtrWV
	int* rowPtrWV_d, *colIndWV_d; dtype* WV_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrWV_d, (num_colJc + 1)*sizeof(rowPtrWV_d[0]));

	if ((cudaStat1 == cudaSuccess)){
		cout << "34.Successfully Allocated CSR rowPtr for WV in GPU memory" << endl;
	}
	else{
		cout << "*** No (34) *** Allocation of CS rowPtr for WV in GPU Memory failed" << endl;
	}

	cudaDeviceSynchronize();
	int jWVId, baseWV;
	int* nnzTotalDevHostPtrWV = &jWVId;
	status2 = cusparseXcsrgemmNnz(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJp, num_colJp, descr, jWId, rowPtrW_d, colIndW_d, descr, jVId, rowPtrV_d, colIndV_d, descr, rowPtrWV_d, nnzTotalDevHostPtrWV);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "35.NNz in WV is done successfully" << endl;
	}
	else{
		cout << "*** No (35) *** NNz in WV failed" << endl;
	}

	cudaDeviceSynchronize();
	if (NULL != nnzTotalDevHostPtrWV){
		jWVId = *nnzTotalDevHostPtrWV;
	}
	else{
		cudaMemcpy(&jWVId, rowPtrWV_d + num_colJc, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseWV, rowPtrWV_d, sizeof(int), cudaMemcpyDeviceToHost);
		jWVId -= baseWV;
	}

	/*cout << "38.The values of nnzWV is " << jWVId << " and total are " << ((15 * num_cam) * (3 * num_pt)) << endl;
	if (jWId == ((15 * num_cam) * (3 * num_pt))){
		cout << "26.*** All W elements are filled" << endl;
	}
	else{
		cout << "26.*** All W are not filled" << endl;
	}*/

	cudaStat1 = cudaMalloc((void**)&colIndWV_d, jWVId*sizeof(colIndWV_d[0]));
	cudaStat2 = cudaMalloc((void**)&WV_d, jWId*sizeof(WV_d[0]));
	if ((cudaStat1 == cudaSuccess) && (cudaStat2 == cudaSuccess)){
		cout << "36.Successfully initialized WV and colIndWV" << endl;
	}
	else{
		cout << "*** No (36) *** Initialization of WV and colIndWV failed" << endl;
	}

	cudaDeviceSynchronize();
	//Finally calculating U
	status2 = cusparseDcsrgemm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJp, num_colJp, descr, jWId, W_d, rowPtrW_d, colIndW_d, descr, jVId, V_d, rowPtrV_d, colIndV_d, descr, WV_d, rowPtrWV_d, colIndWV_d);
	if ((status2 == CUSPARSE_STATUS_SUCCESS)){
		cout << "37.WV is calculated successfully" << endl;
	}
	else{
		cout << "*** No (37) *** WV calculation failed" << endl;
	}

	cudaDeviceSynchronize();
	cout << "************************************************************************" << endl;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Free V
	cudaFree(V_d);
	cudaFree(rowPtrV_d);
	cudaFree(colIndV_d);

	////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate WV * W'
	//Initialize and Allocate rowPtrWVW'
	int* rowPtrWVWT_d, *colIndWVWT_d; dtype* WVWT_d;
	cudaStat1 = cudaMalloc((void**)&rowPtrWVWT_d, (num_colJc + 1)*sizeof(rowPtrWVWT_d[0]));

	if ((cudaStat1 == cudaSuccess)){
		cout << "38.Successfully Allocated CSR rowPtr for WVWT in GPU memory" << endl;
	}
	else{
		cout << "*** No (38) *** CSR rowPtr for WVWT allocation failed" << endl;
	}

	//Create WT in CSR format. WT in CSR == W in CSC
	//Initialize and Allocate csrWT
	dtype* WT_d; int* rowPtrWT_d, *colIndWT_d;
	cudaStat2 = cudaMalloc((void**)&rowPtrWT_d, (num_colJp + 1)*sizeof(rowPtrWT_d[0]));
	cudaStat3 = cudaMalloc((void**)&WT_d, jWId*sizeof(WT_d[0]));
	cudaStat4 = cudaMalloc((void**)&colIndWT_d, jWId*sizeof(colIndWT_d[0]));
	if ((cudaStat2 == cudaSuccess) && (cudaStat3 == cudaSuccess) && (cudaStat4 == cudaSuccess)){
		cout << "39.Successfully initialized WT in GPU memory" << endl;
	}
	else{
		cout << "*** No (39) *** WT initialization in GPU Memory failed" << endl;
	}
	cudaDeviceSynchronize();
	status4 = cusparseDcsr2csc(handle, num_colJc, num_colJp, jWId, W_d, rowPtrW_d, colIndW_d, WT_d, colIndWT_d, rowPtrWT_d, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO);
	if ((status4 == CUSPARSE_STATUS_SUCCESS)){
		cout << "40.CSR to CSC conversion done successfully for WT" << endl;
	}
	else{
		cout << "*** No (40) *** CSR to CSC conversion failed" << endl;
	}

	cudaDeviceSynchronize();
	int jWVWTId, baseWVWT;
	int* nnzTotalDevHostPtrWVWT = &jWVWTId;
	status5 = cusparseXcsrgemmNnz(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJc, num_colJp, descr, jWVId, rowPtrWV_d, colIndWV_d, descr, jWId, rowPtrWT_d, colIndWT_d, descr, rowPtrWVWT_d, nnzTotalDevHostPtrWVWT);
	if ((status5 == CUSPARSE_STATUS_SUCCESS)){
		cout << "41.NNz in WVWT is done successfully" << endl;
	}
	else{
		cout << "*** No (41) *** NNz in WVWT failed" << endl;
	}

	cudaDeviceSynchronize();
	if (NULL != nnzTotalDevHostPtrWVWT){
		jWVWTId = *nnzTotalDevHostPtrWVWT;
	}
	else{
		cudaMemcpy(&jWVWTId, rowPtrWVWT_d + num_colJc, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&baseWVWT, rowPtrWVWT_d, sizeof(int), cudaMemcpyDeviceToHost);
		jWVWTId -= baseWVWT;
	}

	cout << "The values of nnzWVWT is " << jWVWTId << " and total are " << (num_cam * num_cam * 225) << endl;
	/*if (jUId == (num_cam * 225)){
		cout << "12.*** All the U matrices are filled" << endl;
	}
	else{
		cout << "12.*** All U are not filled" << endl;
	}*/

	cudaStat5 = cudaMalloc((void**)&colIndWVWT_d, jWVWTId*sizeof(colIndWVWT_d[0]));
	cudaStat6 = cudaMalloc((void**)&WVWT_d, jWVWTId*sizeof(WVWT_d[0]));
	if ((cudaStat5 == cudaSuccess) && (cudaStat6 == cudaSuccess)){
		cout << "42.Successfully initialized WVWT and colIndWVWT" << endl;
	}
	else{
		cout << "*** No (42) *** Initialization of WVWT and colIndWVWT failed" << endl;
	}

	cudaDeviceSynchronize();
	//Finally calculating U
	status6 = cusparseDcsrgemm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, num_colJc, num_colJc, num_colJp, descr, jWVId, WV_d, rowPtrWV_d, colIndWV_d, descr, jWId, WT_d, rowPtrWT_d, colIndWT_d, descr, WVWT_d, rowPtrWVWT_d, colIndWVWT_d);
	if ((status6 == CUSPARSE_STATUS_SUCCESS)){
		cout << "43.WVWT is calculated successfully" << endl;
	}
	else{
		cout << "*** No (43) *** WVWT calculation failed" << endl;
	}
	cudaDeviceSynchronize();
	cout << "************************************************************************" << endl;
	////////////////////////////////////////////////////////////////////////////////////////////////
	//Free Data
	cudaFree(WV_d);
	cudaFree(rowPtrWV_d);
	cudaFree(colIndWV_d);
	cudaFree(W_d);
	cudaFree(rowPtrW_d);
	cudaFree(colIndW_d);
	cudaFree(WT_d);
	cudaFree(rowPtrWT_d);
	cudaFree(colIndWT_d);
	////////////////////////////////////////////////////////////////////////////////////////////////
	//Calculate S = U - WVWT
	sCalculation(U_d, rowPtrU_d, colIndU_d, WVWT_d, rowPtrWVWT_d, colIndWVWT_d);

auto end = chrono::high_resolution_clock::now();
        cout << "Time taken for S calculation in milliseconds is " << (chrono::duration_cast<chrono::milliseconds>(end - start).count()) << endl;

	////////////////////////////////////////////////////////////////////////////////////////////////
	cudaFree(U_d);
	cudaFree(rowPtrU_d);
	cudaFree(colIndU_d);
	cudaFree(WVWT_d);
	cudaFree(rowPtrWVWT_d);
	cudaFree(colIndWVWT_d);
	////////////////////////////////////////////////////////////////////////////////////////////////
	//Remove matrix descriptor
	status1 = cusparseDestroyMatDescr(descr);
	descr = 0;
	if (status1 == CUSPARSE_STATUS_SUCCESS){
		cout << "44.Descriptor Destroyed" << endl;
	}
	else{
		cout << "*** No (44) *** Destruction of descriptor failed" << endl;
	}

	//Remove handle
	status1 = cusparseDestroy(handle);
	handle = 0;
	if (status1 == CUSPARSE_STATUS_SUCCESS){
		cout << "45.Handle released" << endl;
	}
	else{
		cout << "*** No (45) *** Handle release failed" << endl;
	}
	cout << "************************************************************************" << endl;
	cout << "************************************************************************" << endl;
}

dtype muCalculation()
{
	dtype uMax, vMax;

	uMax = maxUDiagonal();
	vMax = maxVDiagonal();

	if (uMax >= vMax) return (to*uMax);
	return (to*vMax);
}

void augmentUV()
{
	int camId, ptId;
	//Augment U
	for (int i = 0; i < num_cam; i++){
		for (int j = 0; j < 15; j++){
			camId = (i * 225) + (16 * j);
			U[camId] = U[camId] + mu;
		}
	}
	//Augment V
	for (int i = 0; i < num_pt; i++){
		for (int j = 0; j < 3; j++){
			ptId = (i * 9) + (4 * j);
			V[ptId] = V[ptId] + mu;
		}
	}
}
