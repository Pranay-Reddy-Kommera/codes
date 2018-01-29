/////////////////////////////////////////////////////////////////////////////
///   File           :            BA_declaration.h
///   Description    :
/////////////////////////////////////////////////////////////////////////////

#pragma once 

#include "BA_include.h"

/*-------------------------------------------------------------------------*/
/*                      Global variables declaration                       */
/*-------------------------------------------------------------------------*/
extern unsigned int num_cam, num_pt, num_proj;
extern dtype *P, *projections, *calcProjections, *Pnew, *calcProjectionsnew;//, *projectionsFile;
extern unsigned int *ptidx, *camidx;//, *ptidxFile, *camidxFile;
extern dtype *Ep, *Epnew;
extern vector<dtype> errorVector, cgTime, cgIter;
extern dtype *U, *V, *b, *M;
extern dtype *valW, mu;
//extern unsigned int *rowPtrW;
extern bool stop;
extern dtype *delta;
extern dtype Epvector, Epvectornew, rho, vv;

//Constants
const double e1 = pow(10, -25);
const double to = pow(10, -3);
const int kmax = 5;
const double e2 = 0;

/*-------------------------------------------------------------------------*/
/*                           Function declaration                          */
//-------------------------------------------------------------------------//
//Read data from the provided text file
bool readData(char* fileName);

//Declare Dimensions
void declareDim();

//Read complete data from the file
bool read(ifstream& infile);

//Change format
void changeFormat();

//Write input data into text files
bool writeInputData();

//Rodrigues rotation
void setRotationMatrix(int i, dtype *k, int idx);

//Calculate vector norm
dtype vecNorm(dtype *k, int size);

//Clear data
void clearData();

//Calculate measured projections
void measuredProjections();

//Calculate error vector
void errorVectorCalculation();

//Calculate rms value of error
void errorRMSCalculation();

//Calculate A and b
void calculateAb();

//Normal equation
void normalEquations(jacobian* J);

//Initialize U, V and W with zeros
void zeroInitialization();

//Calculate U,V and W
void calculateUVW(jacobian* J);

//Calculate U vector
void matMulU(dtype Aij[2][15], int cam);

//Calculate V vector
void matMulV(dtype Bij[2][3], int pt);

//Calculate W vector
void matMulW(dtype Aij[2][15], dtype Bij[2][3], int proj);

//Calculate Ea vector
void matMulEa(dtype Aij[2][15], int i);

//Calculate Eb vector
void matMulEb(dtype Bij[2][3], int i);

//Calculate b
void calculateb(jacobian* J);

//Calculate stoppign criteria
bool stopCriteria();
bool criteria();

//Damping term calculation
dtype muCalculation();
dtype maxUDiagonal();
dtype maxVDiagonal();

//Augment U and V
void augmentUV();

//Delta initialization
void deltaInitialization();

//Conjugate gradient without preconditioner
void cgwjp();

//rk = b
void rkEqualsb(dtype* rk);

//Compute zk = M^-1 * rk
void computezk(dtype* zk, dtype* rk);

//pk = zk
void pkEqualszk(dtype* pk, dtype* zk);

//sk' * sk
dtype vecMulVec(dtype* sk);
dtype vecMulVec(dtype* pk, dtype* Apk);

//A * pk
void matMulVec(dtype* Apk, dtype* pk);

//Apk initialization
void apkInitialization(dtype* Apk);

//delta = delta + (alphak * pk)
void updatedelta(dtype alphak, dtype* pk);

//rk = rk - (alphak * Apk)
void updaterk(dtype* rk, dtype alphak, dtype* Apk);
void updaterkp1(dtype*rkp1, dtype* rk, dtype alphak, dtype* Apk);

//pk = zk + (betak * pk)
void updatepk(dtype* pk, dtype* rk, dtype* zk, dtype zktrk);
void updatepk(dtype* pk, dtype* rkp1, dtype* rk, dtype* zk, dtype zktrk);

//rk = rkp1
void copyrkp1Tork(dtype* rkp1, dtype* rk);

//Print CG details
void printCGDetails();

//Update parameters Pnew
void updateParametersPnew();

//Calculate new f(Pnew)
void measuredProjectionsnew();

//Calculate new error vector
void EpnewCalculation();

//Calculate denominator
dtype calculateDenom();

//Update parameter P
void updateParametersP();

//Print error vector
void printErrorVector();

//Compute Jacobian Preconditioner
void computeJacobiPreconditioner();

//-------------------------------------------------------------------------//
