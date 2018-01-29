///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_includeFiles.h
///   Description    :         Includes all necessary library files and declarations.
//////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <ctime>
#include <chrono>
#include "cholmod.h"
#include <stdlib.h>

using namespace std;

typedef double dtype;

//**********************************************************
//Data types and Data declarations.
//typedef double dtype;

extern int num_cam, num_pt, num_proj;
extern dtype **R, **T, *focal, **rd, *points, *projections, *P, *calcProjections, *Pnew, *calcProjectionsnew, *Ep, *Epnew;
extern int *ptidx, *camidx, *ncpoints;
extern dtype *U, *V, *W, *g, mu, *delta, *Vinv;
extern dtype rho, Epvector, Epvectornew, vv;

extern vector<dtype> errorVector, rhoVector, muVector;
extern bool stop;

//Constants
const double e1 = pow(10, -12);
const double e3 = pow(10, -50);
const double e2 = 0;
const double to = pow(10, -3);
const int kmax = 5;

//A and B declarations
struct Jacobian{
	dtype Aij[2][15];
	dtype Bij[2][3];
};

//Intermediate S values
struct Sint{
	dtype Sij[15][15];
};

//**********************************************************

