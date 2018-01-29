/////////////////////////////////////////////////////////////////////////////
///   File           :            BA_include.h
///   Description    :
/////////////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>
#include <vector>
#include <iomanip>

using namespace std;

typedef double dtype;

//Jacobian structure
struct jacobian{
   dtype Aij[2][15];
   dtype Bij[2][3];
};

//#define _write
//#define Polak_Ribiere
