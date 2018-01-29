///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_Math.cpp
///   Description    :         Simple mathematical operations.
//////////////////////////////////////////////////////////////////////////////////////

#include "BA_includeFiles.h"
#include "BA_data.h"

/**********************************************************************/
/*                Calculate L2 vector Norm                            */
/**********************************************************************/
dtype vecNorm(dtype *k, int size){
   dtype tmp = 0;

   for (int j = 0; j < size; j++)
   {
      tmp += k[j] * k[j];
   }

   return (sqrt(tmp));
}

/**********************************************************************/
/*                    Calculate U                                     */
/**********************************************************************/
void MatMulU(dtype Aij[2][15], int cam)
{
   dtype x, y; int id;

   for (int i = 0; i < 15; i++){
      x = Aij[0][i]; y = Aij[1][i];
      for (int j = 0; j < 15; j++){
         id = (i * 15) + j + cam;
         U[id] = U[id] + (x * Aij[0][j]) + (y * Aij[1][j]);
      }
   }
}

/**********************************************************************/
/*                             Calculate V                            */
/**********************************************************************/
void MatMulV(dtype Bij[2][3], int pt)
{
   dtype x, y; int id;

   for (int i = 0; i < 3; i++){
      x = Bij[0][i]; y = Bij[1][i];
      for (int j = 0; j < 3; j++){
         id = (i * 3) + j + pt;
         V[id] = V[id] + (x * Bij[0][j]) + (y * Bij[1][j]);
      }
   }
}

/**********************************************************************/
/*                   Calculate W                                      */
/**********************************************************************/
void MatMulW(dtype Aij[2][15], dtype Bij[2][3], int proj)
{
   dtype x, y;

   for (int i = 0; i < 15; i++){
      x = Aij[0][i]; y = Aij[1][i];
      for (int j = 0; j < 3; j++){
         W[(i * 3) + j + proj] = (x * Bij[0][j]) + (y * Bij[1][j]);
      }
   }
}

/**********************************************************************/
/*                      Calculate Ea                                  */
/**********************************************************************/
void MatMulEa(dtype Aij[2][15], int i)
{
   int camId, Id; dtype x, y;
   camId = (camidx[i] * 15); Id = 2 * i;
   //x = projections[Id] - calcProjections[Id];
   //y = projections[Id + 1] - calcProjections[Id + 1];

   for (int j = 0; j < 15; j++){
      //g[camId + j] = g[camId + j] + (x*Aij[0][j]) + (y*Aij[1][j]);
      g[camId + j] = g[camId + j] + (Ep[Id]*Aij[0][j]) + (Ep[Id+1]*Aij[1][j]);
   }
}

/**********************************************************************/
/*                             Calculate Eb                           */
/**********************************************************************/
void MatMulEb(dtype Bij[2][3], int i)
{
   int ptId, Id; dtype x, y;
   ptId = (num_cam * 15) + (ptidx[i] * 3); Id = 2 * i;
   //x = projections[Id] - calcProjections[Id];
   //y = projections[Id + 1] - calcProjections[Id + 1];

   for (int j = 0; j < 3; j++){
      //g[ptId + j] = g[ptId + j] + (x*Bij[0][j]) + (y*Bij[1][j]);
      g[ptId + j] = g[ptId + j] + (Ep[Id]*Bij[0][j]) + (Ep[Id+1]*Bij[1][j]);
   }
}

/*void invCalculation(dtype V[], dtype* Vinv)
{
   dtype a, b, c, d, e, f, g, h, i, det;
   a = V[0]; b = V[1]; c = V[2];
   d = V[3]; e = V[4]; f = V[5];
   g = V[6]; h = V[7]; i = V[8];

   det = (a*((e*i) - (h*f))) - (b*((d*i) - (g*f))) + (c*((d*h) - (e*g)));

   Vinv[0] = ((e*i) - (h*f)) / det; Vinv[1] = ((h*c) - (b*i)) / det; Vinv[2] = ((b*f) - (e*c)) / det;
   Vinv[3] = ((f*g) - (d*i)) / det; Vinv[4] = ((a*i) - (g*c)) / det; Vinv[5] = ((d*c) - (a*f)) / det;
   Vinv[6] = ((d*h) - (e*g)) / det; Vinv[7] = ((b*g) - (a*h)) / det; Vinv[8] = ((a*e) - (b*d)) / det;

}*/

/**********************************************************************/
/*                        Calculate Vinv                              */
/**********************************************************************/
void invCalculation(dtype* V, dtype* Vinv)
{
   dtype a, b, c, d, e, f, g, h, i, det;
   a = V[0]; b = V[1]; c = V[2];
   d = V[3]; e = V[4]; f = V[5];
   g = V[6]; h = V[7]; i = V[8];

   det = (a*((e*i) - (h*f))) - (b*((d*i) - (g*f))) + (c*((d*h) - (e*g)));

   Vinv[0] = ((e*i) - (h*f)) / det; Vinv[1] = ((h*c) - (b*i)) / det; Vinv[2] = ((b*f) - (e*c)) / det;
   Vinv[3] = ((f*g) - (d*i)) / det; Vinv[4] = ((a*i) - (g*c)) / det; Vinv[5] = ((d*c) - (a*f)) / det;
   Vinv[6] = ((d*h) - (e*g)) / det; Vinv[7] = ((b*g) - (a*h)) / det; Vinv[8] = ((a*e) - (b*d)) / det;
}

/**********************************************************************/
/*                        Calculate Y                                 */
/**********************************************************************/
void MatMulY(dtype* Yij, int i)
{
   dtype tmp;

   for (int j = 0; j < 15; j++){
      for (int k = 0; k < 3; k++){
         tmp = 0;
         for (int l = 0; l < 3; l++){
            tmp = tmp + (W[(i * 45) + (j * 3) + l] * Vinv[(ptidx[i] * 9) + k + (l * 3)]);
         }
         Yij[(i * 45) + (j * 3) + k] = tmp;
      }
   }
}

/**********************************************************************/
/*                 Calculate intermediate S values                    */
/**********************************************************************/
void interSCalculation(dtype Yij[], int j, Sint* stmp)
{
   dtype tmp;

   for (int i = 0; i < 15; i++){
      for (int k = 0; k < 15; k++){
         tmp = 0;
         for (int l = 0; l < 3; l++){
            tmp = tmp + (Yij[(i * 3) + l] * W[(j * 45) + (k * 3) + l]);
         }
         stmp->Sij[i][k] = stmp->Sij[i][k] + tmp;
      }
   }

}

/**********************************************************************/
/*                 Calculate intermediate rhs value                   */
/**********************************************************************/
void interhsCalculation(dtype* Yij, dtype* rhs, int camId, int ptId)
{
   dtype tmp; int Id = num_cam * 15;

   for (int i = 0; i < 15; i++){
      tmp = 0;
      for (int j = 0; j < 3; j++){
         tmp = tmp + (Yij[(i * 3) + j] * g[Id + (ptId * 3) + j]);
      }
      rhs[(camId * 15) + i] = rhs[(camId * 15) + i] + tmp;
   }
}

/**********************************************************************/
/*                   Cholesky Decomposition                           */
/**********************************************************************/
void chol_decom(dtype** S, dtype* l, int n)
{
   dtype tmp;

   for (int i = 0; i < n; i++){
      for (int j = 0; j <= i; j++){
         tmp = 0;
         for (int k = 0; k < j; k++){
            tmp += l[((i*(i + 1)) / 2) + k] * l[((j*(j + 1)) / 2) + k];
         }
         if (i == j){
            l[((i*(i + 1)) / 2) + j] = sqrt(S[i][j] - tmp);
         }
         else{
            l[((i*(i + 1)) / 2) + j] = (1 / l[((j*(j + 3)) / 2)])*(S[i][j] - tmp);
         }
      }
   }
   //cout << l[49760] << "\t" << l[49767] << "\t" << l[49769] << endl;
}

/**********************************************************************/
/*                  Forward calculation                               */
/**********************************************************************/
void forwardCalculation(dtype* l, dtype* y, dtype* rhs)
{
   dtype tmp;

   for (int i = 0; i < (num_cam * 15); i++){
      tmp = 0;
      for (int j = 0; j < i; j++){
         tmp += (l[((i*(i + 1)) / 2) + j] * y[j]);
      }
      y[i] = (1 / l[((i*(i + 3)) / 2)])*(rhs[i] - tmp);
   }
}

/**********************************************************************/
/*                         Backward Calculation                       */
/**********************************************************************/
void backwardCalculation(dtype* l, dtype* y)
{
   dtype tmp;

   for (int i = (num_cam * 15) - 1; i >= 0; i--){
      tmp = 0;
      for (int j = i + 1; j < (num_cam * 15); j++){
         tmp += (l[((j*(j + 1)) / 2) + i] * delta[j]);
      }
      delta[i] = (1 / l[((i*(i + 3)) / 2)])*(y[i] - tmp);
   }

   //cout << delta[63] << "\t" << delta[157] << "\t" << delta[249] << endl;
}

/**********************************************************************/
/*                         Calculate deltaB                           */
/**********************************************************************/
void MatMuldeltaB(int i)
{
   int camId, ptId, projId, Id;
   camId = camidx[i]; ptId = ptidx[i]; projId = i * 45; Id = num_cam * 15;
   dtype tmp;
	
   for (int j = 0; j < 3; j++){
      tmp = 0;
      for (int k = 0; k < 15; k++){
         tmp += (W[projId + (k * 3) + j])*(delta[(camId * 15) + k]);
      }
      delta[Id + (ptId * 3) + j] += tmp;
   }

}

/**********************************************************************/
/*                       Calculate deltaB final                       */
/**********************************************************************/
void MatMuldeltaBFinal(dtype* Vinv, int i)
{
   dtype tmp[3];
   int Id = num_cam * 15; int Idd = Id + (i * 3);

   tmp[0] = (Vinv[0] * delta[Idd]) + (Vinv[1] * delta[Idd + 1]) + (Vinv[2] * delta[Idd + 2]);
   tmp[1] = (Vinv[3] * delta[Idd]) + (Vinv[4] * delta[Idd + 1]) + (Vinv[5] * delta[Idd + 2]);
   tmp[2] = (Vinv[6] * delta[Idd]) + (Vinv[7] * delta[Idd + 1]) + (Vinv[8] * delta[Idd + 2]);

   delta[Idd] = tmp[0];
   delta[Idd + 1] = tmp[1];
   delta[Idd + 2] = tmp[2];
}

/**********************************************************************/
