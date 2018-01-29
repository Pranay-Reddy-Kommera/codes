///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_functions.cpp
///   Description    :         Functions required in Bundle Adjustment
//////////////////////////////////////////////////////////////////////////////////////

#include "BA_includeFiles.h"
#include "BA_data.h"

/*****************************************************************/
/*               Calculate Measured projections                  */
/*****************************************************************/
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

/*****************************************************************/
/*        Calculate Measured Projections for new parameters      */
/*****************************************************************/
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

/*****************************************************************/
/*                     Calculate error vector                    */
/*****************************************************************/
void errorCalculation()
{
   for (int i = 0; i < (2 * num_proj); i++){
      Ep[i] = projections[i] - calcProjections[i];
   }
}

/*****************************************************************/
/*                     Calculate RMS error Values                */
/*****************************************************************/
void RMSErrorCalculation()
{
   dtype EpNorm, tmp, rms;

   EpNorm = vecNorm(Ep, 2*num_proj);
   tmp = (EpNorm*EpNorm) / (2 * num_proj);
   rms = sqrt(tmp);
   errorVector.push_back(rms);

   //cout << "RMS Error Norm is " << rms << endl;

}

/*****************************************************************/
/*                   Calculate U, V, W, Ea, Eb                   */
/*****************************************************************/
void UVWgCalculation()
{
   //Allocate Memory for Jacobian
   Jacobian* J = new Jacobian[num_proj];

   //Calculate Normal equation parameters
   NormalEquations(J);

   //cout << J[10].Aij[0][0] << "\t" << J[10].Bij[1][1] << endl;

   //Initialize U, V, W to zeros
   UVWgZeroInit();

   //cout << U[450] << "\t" << V[18] << "\t" << W[90] << endl;
   //cout << g[15] << "\t" << g[16] << "\t" << g[318] << "\t" << g[319] << endl;

   //Calculate U, V, W
   UVW(J);

   //cout << U[450] << "\t" << V[18] << "\t" << W[90] << endl;

   //Calculate g=[Ea',Eb']
   EaEb(J);

   //cout << g[15] << "\t" << g[16] << "\t" << g[318] << "\t" << g[319] << endl;

   //Delete Jacobian
   delete[] J;
}

/*****************************************************************/
/*                        Stop Criteria                          */
/*****************************************************************/
bool stopCriteria()
{
   return criteria();
}

/*****************************************************************/
/*                    Calculate damping term                     */
/*****************************************************************/
dtype muCalculation()
{
   dtype uMax, vMax;

   uMax = maxUDiagonal();
   vMax = maxVDiagonal();

   if (uMax >= vMax) return (to*uMax);
   return (to*vMax);
}

/*****************************************************************/
/*                    Delta Calculation                          */
/*****************************************************************/
void deltaCalculation()
{
   //int* ncpoints = new int[num_pt];
	
   //S = [num_cam*15][num_cam*15]
   //Memory for row elements
   dtype** S = new dtype *[num_cam * 15];
   //Memory for column elements
   for (int i = 0; i < (num_cam * 15); i++){
      S[i] = new dtype[num_cam * 15];
   }

   //rhs = [nu_cam*15]
   dtype* rhs = new dtype[num_cam * 15];

   //Calculate augmented U and V matrix
   augmentUV();

   //cout<<"Augmentation is done..."<<endl;

   //Create Y matrix
   dtype* Yij = new dtype[num_proj * 45];

   //auto start1 = chrono::high_resolution_clock::now();
   //Calculate Yij matrix
   invVCalculation();
   //cout<<"inv is calculated"<<endl;
   yCalculation(Yij);

   //cout<<"Yij done..."<<endl;

   //Calculate the number of camers each points are in
   //numCamPtsCalculations(ncpoints);

   //cout << ncpoints[11314] << "\t" << ncpoints[11000] << "\t" << ncpoints[6000] << endl;

   //Compute Matrix S
   sCalculation(S, Yij, ncpoints);

   //auto end1 = chrono::high_resolution_clock::now();
   //cout << "Time taken for S in milliseconds is " << chrono::duration_cast<chrono::milliseconds>(end1 - start1).count() << endl;
   //Check the sparsity of S matrix
   int nnz = sparsityCheck(S, num_cam);

   //Read S into ST format
   stAFormat(S, num_cam, nnz);

   //Compute rhs
   rhsCalculation(rhs, Yij);
   //cout<<"rhs is calculated"<<endl;

   //Read file in ST format
   stbFormat(rhs, num_cam);

   //auto start2 = chrono::high_resolution_clock::now();
   //Calculate deltaA using cholesky decomposition
   //deltaACalculation(S, rhs);
   deltaACalculationCHOLMOD("ST_AFormat.txt","ST_bFormat.txt","ST_xFormat.txt");
   //auto end2 = chrono::high_resolution_clock::now();
   //cout << "Time taken for cholmod  in milliseconds is " << chrono::duration_cast<chrono::milliseconds>(end2 - start2).count() << endl;
   //Retrieve x from CHOLMD
   extractX("ST_xFormat.txt");

   //cout<<"Chol_mod is calculated"<<endl;

   //Calculate deltaB
   deltaBCalculation();

   //Delete local variables
   //S = [num_cam*15][num_cam*15]
   for (int i = 0; i < (num_cam*15); i++){
      delete[] S[i];
   }
   delete[] S;
   delete[] rhs;
   //delete[] ncpoints;
   delete[] Yij;
}

/*****************************************************************/
/*                    Update parameter Pnew                      */
/*****************************************************************/
void updateParametersPnew()
{
   for (int i = 0; i < ((num_cam * 15) + (num_pt * 3)); i++){
      Pnew[i] = P[i] + delta[i];
   }
}

/*****************************************************************/
/*                  New error vector calculation                 */
/*****************************************************************/
void EpnewCalculation()
{
   for (int i = 0; i < 2 * num_proj; i++){
      Epnew[i] = projections[i] - calcProjectionsnew[i];
   }
}

/*****************************************************************/
/*                   Calculate denominator                       */
/*****************************************************************/
dtype calculateDenom()
{
   int Id = (num_cam * 15) + (num_pt * 3);
   dtype *tmp = new dtype[Id];
   dtype denom = 0;

   for (int i = 0; i < Id; i++){
      tmp[i] = (mu*delta[i]) + g[i];
   }

   for (int i = 0; i < Id; i++){
      denom += delta[i] * tmp[i];
   }
   //cout << "denom = " << denom << endl;
   return denom;
}

/*****************************************************************/
/*                    Update parameter P                         */
/*****************************************************************/
void updateParametersP()
{
   for (int i = 0; i < ((num_cam * 15) + (num_pt * 3)); i++){
      P[i] = Pnew[i];
   }
}

/*****************************************************************/
/*                     Update CalcProjections                    */
/*****************************************************************/
void updateCalcProjections()
{
   for (int i = 0; i < (2 * num_proj); i++){
      calcProjections[i] = calcProjectionsnew[i];
   }
}

/*****************************************************************/
/*                      Update Ep                                */
/*****************************************************************/
void updateEp()
{
   for (int i = 0; i < (2 * num_proj); i++){
      Ep[i] = Epnew[i];
   }
}

/*****************************************************************/
