/////////////////////////////////////////////////////////////////////////////
///   File           :          BA_functions.cpp
///   Description    :
/////////////////////////////////////////////////////////////////////////////

#include "BA_declaration.h"

/*-------------------------------------------------------------------------*/
/*      Change camidx, ptidx and projections into cam format               */
/*-------------------------------------------------------------------------*/
/*void changeFormat()
{
   int Id = 0;

   for(int cam=0;cam<num_cam;cam++){
      for(int i=0;i<num_proj;i++){
         if(camidxFile[i]==cam){
            int Idx = Id*2;
            camidx[Id]=cam; ptidx[Id]=ptidxFile[i];
            projections[Idx]=projectionsFile[i*2]; 
            projections[Idx+1]=projectionsFile[(i*2)+1];
            ++Id;
         }
      }
   }

}
*/
/*-------------------------------------------------------------------------*/
/*                      Calculate Measured projections                     */
/*-------------------------------------------------------------------------*/
//Calculating f(P)
void measuredProjections()
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

/*-------------------------------------------------------------------------*/
/*                        Calculate error vector                           */
/*-------------------------------------------------------------------------*/
void errorVectorCalculation()
{
   for(int i=0;i<(2*num_proj);i++){
      Ep[i]=projections[i]-calcProjections[i];
   }
}

/*-------------------------------------------------------------------------*/
/*                   Calculate rms values of error vector                  */
/*-------------------------------------------------------------------------*/
void errorRMSCalculation()
{
   dtype EpNorm, tmp, rms;

   EpNorm = vecNorm(Ep, 2*num_proj);
   tmp = (EpNorm*EpNorm) / (2*num_proj);
   rms = sqrt(tmp);
   
   errorVector.push_back(rms);
//   cout<<"RMS error is "<<rms<<endl;
}

/*-------------------------------------------------------------------------*/
/*                           Calculate A and b                             */
/*-------------------------------------------------------------------------*/
void calculateAb()
{
   long long unsigned int totalRows, totalCols, nzRows;

   //Allocate memory for Jacobian
   jacobian* J = new jacobian[num_proj];

   //Calculate normal equation
   normalEquations(J);
   totalRows = num_pt*num_cam*2;
   totalCols = (15*num_cam)+(3*num_pt);
   nzRows = num_proj*2;
//   cout<<"J of size "<<totalRows<<" x "<<totalCols<<" with "<<nzRows<<" non-zero rows is calculated"<<endl;

   //initializeUVW with zeros
   zeroInitialization();

   //Calculate UVW
   calculateUVW(J);

   //Calculate b=[Ea' Eb']
   calculateb(J);

   //Delete Jacobian
   delete[] J;
}

/*-------------------------------------------------------------------------*/
/*                           Calculate stopping criteria                   */
/*-------------------------------------------------------------------------*/
bool stopCriteria()
{
   return criteria();
}

/*-------------------------------------------------------------------------*/
/*                          Calculate damping term                         */
/*-------------------------------------------------------------------------*/
dtype muCalculation()
{
   dtype uMax, vMax;

   uMax = maxUDiagonal();
   vMax = maxVDiagonal();

   if (uMax >= vMax) return (to*uMax);
   return (to*vMax);
}

/*-------------------------------------------------------------------------*/
/*                       Augment U and V                                   */
/*-------------------------------------------------------------------------*/
void augmentUV()
{
   int camId,ptId;
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

/*-------------------------------------------------------------------------*/
/*              Conjugate gradient without preconditioner                  */
/*-------------------------------------------------------------------------*/
void cgwjp()
{
   //Initialize delta to 0; xk = 0
   deltaInitialization();

   int size = (num_cam*15)+(num_pt*3);
   //Initialize sk, pk ,Apk, sktsk and k
   dtype* rk = new dtype[size];
#if defined Polak_Ribiere
   dtype* rkp1 = new dtype[size];
#endif
   dtype* zk = new dtype[size];
   dtype* pk = new dtype[size];
   dtype* Apk = new dtype[size];
   dtype zktrk, alphak;
   int kk;
   
   //rk = b - (A * xk);as xk = 0 rk = b
   rkEqualsb(rk);
   //zk = M^-1 * rk
   computezk(zk,rk);
   //pk = szk
   pkEqualszk(pk,zk);
   kk = 0;

   //Start time stamp
   auto start1 = chrono::high_resolution_clock::now();

   //While loop
   while((vecNorm(rk,size)/vecNorm(b,size))>=0.001){

      apkInitialization(Apk);
      //zk' * rk
      zktrk = vecMulVec(zk,rk);
      //Apk = A * pk
      matMulVec(Apk, pk);
      //alphak = zktrk/(pk' * Apk)
      alphak = zktrk/(vecMulVec(pk,Apk));
      //delta = delta + alphak * pk
      updatedelta(alphak, pk);
#if defined Polak_Ribiere
      //rkp1 = rk - (alphak * Apk)
      updaterkp1(rkp1, rk, alphak, Apk);
      //zk = M^-1 * rkp1
      computezk(zk,rkp1);
      //pk = rk + (betak * pk)
      updatepk(pk,rkp1,rk,zk,zktrk);
#else
      //rk = rk - (alphak * Apk)
      updaterk(rk, alphak, Apk);
      //zk = M^-1 * rk
      computezk(zk,rk);
      //pk = rk + (betak * pk)
      updatepk(pk,rk,zk,zktrk);
#endif
      //k = k + 1
      kk = kk + 1;
#if defined Polak_Ribiere
      //rk = rkp1
      copyrkp1Tork(rkp1,rk);
#endif
   }   

   //Stop time
   auto end1 = chrono::high_resolution_clock::now();

   cgIter.push_back(kk);
   cgTime.push_back(chrono::duration_cast<chrono::milliseconds>(end1-start1).count());

   //deallocate sk and pk
   delete[] rk;
#if defined Polak_Ribiere
   delete[] rkp1;
#endif
   delete[] zk;
   delete[] pk;
   delete[] Apk;
}

/*-------------------------------------------------------------------------*/
/*                    Print Conjugate Gradient details                     */
/*-------------------------------------------------------------------------*/
void printCGDetails()
{
   cout<<"\nConjugate Gradient time and iterations"<<endl;
   
   for(int i=0;i<cgTime.size();i++){
      cout<< "Step "<<(i+1)<<": Iterations = "<<cgIter[i]<<" and Time = "<<cgTime[i]<<" milliseconds"<<endl;
   }
}

/*-------------------------------------------------------------------------*/
/*                      Update parameters Pnew                             */
/*-------------------------------------------------------------------------*/
void updateParametersPnew()
{
   int totalSize = (num_cam*15)+(num_pt*3);
   for (int i = 0; i < totalSize; i++){
      Pnew[i] = P[i] + delta[i];
   }
}

/*-------------------------------------------------------------------------*/
/*                  Calculate new f(Pnew)                                  */
/*-------------------------------------------------------------------------*/
void measuredProjectionsnew()
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

/*-------------------------------------------------------------------------*/
/*                    Calculate new error vector Epnew                     */
/*-------------------------------------------------------------------------*/
void EpnewCalculation()
{
   for (int i = 0; i < 2 * num_proj; i++){
      Epnew[i] = projections[i] - calcProjectionsnew[i];
   }
}

/*-------------------------------------------------------------------------*/
/*                  Calculate denominator                                  */
/*-------------------------------------------------------------------------*/
dtype calculateDenom()
{
   int Id = (num_cam * 15) + (num_pt * 3);
   dtype *tmp = new dtype[Id];
   dtype denom = 0;

   for (int i = 0; i < Id; i++){
      tmp[i] = (mu*delta[i]) + b[i];
   }

   for (int i = 0; i < Id; i++){
      denom += delta[i] * tmp[i];
   }
	//cout << "denom = " << denom << endl;
   return denom;
}

/*-------------------------------------------------------------------------*/
/*                     Update parameters P                                 */
/*-------------------------------------------------------------------------*/
void updateParametersP()
{
   int totalSize = (num_cam*15)+(num_pt*3);

   for (int i = 0; i < totalSize; i++){
      P[i] = Pnew[i];
   }
}

/*-------------------------------------------------------------------------*/
/*                          Print error vector                             */
/*-------------------------------------------------------------------------*/
void printErrorVector()
{
   cout << "\nThe error vector is " << endl;

cout<<setprecision(10)<<fixed;

   for (int l = 0; l < errorVector.size(); l++){
      cout << errorVector[l] << endl;
   }
}

/*-------------------------------------------------------------------------*/
/*                     Compute Jacobi Preconditioner                       */
/*-------------------------------------------------------------------------*/
void computeJacobiPreconditioner()
{
   // M is inverse of diag of A=[U V]
   for(int i=0;i<num_cam;i++){
      for(int j=0;j<15;j++){
         M[(i*15)+j]=1/(U[(i*225)+(16*j)]);
      }
   }

   for(int i=0;i<num_pt;i++){
      int Id = (num_cam*15)+(i*3);
      for(int j=0;j<3;j++){
         M[Id+j]=1/(V[(i*9)+(4*j)]);
      }
   }
}

/*-------------------------------------------------------------------------*/
