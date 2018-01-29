///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_data.h
///   Description    :         Files which reads data into respective variables.
//////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "BA_includeFiles.h"

//**********************************************************
//Declarations

//#define _Write Write_Data
#define _PVector Uses_P
//#define _RTVector Uses_RTfrd

//Read data from the text file
bool ReadData(char *fileName);

//Write data to a text files
bool WriteInputData();

//Calculate measured projections
//Calculate f(P)
void MeasuredProjections();
void MeasuredProjectionsnew();

//Calculate RMS Error
//Calculate Ep
void errorCalculation();
void RMSErrorCalculation();

//Calculate JtJ and JtEp
void UVWgCalculation();
void UVWgZeroInit();

//Clear Data
void ClearData();

//Stopping Criteria
bool stopCriteria();

//mu Calculation
dtype muCalculation();

//Rodrigues rotation
void SetRotationMatrix(int i, dtype *k, int idx);

//Calculate vector norm
dtype vecNorm(dtype *k, int size);

//Calculate Jacobian
void NormalEquations(dtype* Jc, dtype* Jp, int* rowJc, int* colJc, int* rowJp, int* colJp);

//U, V, W calculation
void UVW(Jacobian* J);

//Ea, Eb Calculation
void EaEb(Jacobian* J);

//Matrix Multiplication
void MatMulU(dtype Aij[2][15], int cam);
void MatMulV(dtype Bij[2][3], int pt);
void MatMulW(dtype Aij[2][15], dtype Bij[2][3], int proj);
void MatMulEa(dtype Aij[2][15], int i);
void MatMulEb(dtype Bij[2][3], int i);

//Norm Infinity criteria
bool criteria();

//Maximum Diagonal element in U nd V
dtype maxUDiagonal();
dtype maxVDiagonal();

//Calculate delta
void deltaCalculation();

//Augment U and V
void augmentUV();

//Yij Calculation
void invVCalculation();
void inverseVCalculation(dtype* V_d);
void yCalculation(dtype* Yij);

//Inverse
//void invCalculation(dtype V[], dtype* Vinv);
void invCalculation(dtype* V);

//Calculate Yij
void MatMulY(dtype* Yij, int i);

//Calculates the number of camers points are in
//void numCamPtsCalculations(int *ncpoints);
void numCamPtsCalculations();

//S calculation
//void sCalculation();
void sCalculation(dtype* U_d, int* rowPtrU_d, int* colIndU_d, dtype* WVWT_d, int* rowPtrWVWT_d, int* colIndWVWT_d);

//Intermediate S values
void interSCalculation(dtype Yij[], int j, Sint* stmp);

//Check sparsity of matrix
int sparsityCheck(dtype** S, int num_cam);

//Store S in ST format
void stAFormat(dtype** S, int num_cam, int nnz);

//rhs calculation
void rhsCalculation(dtype* rhs, dtype* Yij);

//Store b in ST format
void stbFormat(dtype* rhs, int num_cam);

//Intermediate rhs Calculation
void interhsCalculation(dtype* Yij, dtype* rhs, int camId, int pdId);

//DeltaA Calculation
void deltaACalculation(dtype** S, dtype* rhs);

//Cholesky decomposition
void chol_decom(dtype** S, dtype* l, int n);

//Ly=b calculation
void forwardCalculation(dtype* l, dtype* y, dtype* rhs);

//L'x=y calculation
void backwardCalculation(dtype* l, dtype* y);

//deltaB Calculation
void deltaBCalculation();

//W'deltaA Calculation
void MatMuldeltaB(int i);

//Final deltaB Calculation
void MatMuldeltaBFinal(dtype* Vinv, int i);

//Update parameters Pnew with P
void updateParametersPnew();

//Calculate Ep new
void EpnewCalculation();

//Calculate denominator of rho
dtype calculateDenom();

//Update parameters P with Pnew
void updateParametersP();

//Update calProjections with calcProjectionsnew
void updateCalcProjections();

//Update Ep with Epnew
void updateEp();