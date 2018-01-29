/////////////////////////////////////////////////////////////////////////////
///   File           :          BA_math.cpp
///   Description    :
/////////////////////////////////////////////////////////////////////////////

#include "BA_declaration.h"

/*-------------------------------------------------------------------------*/
/*                  Create Rodrigues rotation from K                       */
/*-------------------------------------------------------------------------*/
void setRotationMatrix(int i, dtype *k, int idx)
{
   //Calculate norm of K
   dtype rNorm, cost, sint;
   rNorm = vecNorm(k, 3);

   if (rNorm == 0){
      cost = 0.5; sint = 1;
   }
   else{
      cost = (1 - cos(rNorm)) / (rNorm*rNorm);
      sint = sin(rNorm) / rNorm;
   }

   P[idx] = 1 - (cost*((k[1] * k[1]) + (k[2] * k[2])));
   P[idx + 1] = (-sint*k[2]) + (cost*k[0] * k[1]);
   P[idx + 2] = (sint*k[1]) + (cost*k[0] * k[2]);
   P[idx + 3] = (sint*k[2]) + (cost*k[0] * k[1]);
   P[idx + 4] = 1 - (cost*((k[0] * k[0]) + (k[2] * k[2])));
   P[idx + 5] = (-sint*k[0]) + (cost*k[1] * k[2]);
   P[idx + 6] = (-sint*k[1]) + (cost*k[0] * k[2]);
   P[idx + 7] = (sint*k[0]) + (cost*k[1] * k[2]);
   P[idx + 8] = 1 - (cost*((k[0] * k[0]) + (k[1] * k[1])));

}

/*-------------------------------------------------------------------------*/
/*                    Calculate L2 norm of a vector k                      */
/*-------------------------------------------------------------------------*/
dtype vecNorm(dtype *k, int size)
{
   dtype tmp = 0;

   for (int j = 0; j < size; j++)
   {
      tmp += k[j] * k[j];
   }

   return (sqrt(tmp));
}

/*-------------------------------------------------------------------------*/
/*                           Calculate vector U                            */
/*-------------------------------------------------------------------------*/
void matMulU(dtype Aij[2][15], int cam)
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

/*-------------------------------------------------------------------------*/
/*                        Calculate vector V                               */
/*-------------------------------------------------------------------------*/
void matMulV(dtype Bij[2][3], int pt)
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

/*-------------------------------------------------------------------------*/
/*                   Calculate vector W                                    */
/*-------------------------------------------------------------------------*/
void matMulW(dtype Aij[2][15], dtype Bij[2][3], int proj)
{
   dtype x, y;

   for (int i = 0; i < 15; i++){
      x = Aij[0][i]; y = Aij[1][i];
      for (int j = 0; j < 3; j++){
         valW[(i * 3) + j + proj] = (x * Bij[0][j]) + (y * Bij[1][j]);
      }
   }
}

/*-------------------------------------------------------------------------*/
/*                         Calculate vector Ea                             */
/*-------------------------------------------------------------------------*/
void matMulEa(dtype Aij[2][15], int i)
{
   int camId, Id; dtype x, y;
   camId = (camidx[i] * 15); Id = 2 * i;

   for (int j = 0; j < 15; j++){
      b[camId + j] = b[camId + j] + (Ep[Id] * Aij[0][j]) + (Ep[Id + 1] * Aij[1][j]);
   }
}

/*-------------------------------------------------------------------------*/
/*                         Calculate vector Eb                             */
/*-------------------------------------------------------------------------*/
void matMulEb(dtype Bij[2][3], int i)
{
   int ptId, Id; dtype x, y;
   ptId = (num_cam * 15) + (ptidx[i] * 3); Id = 2 * i;

   for (int j = 0; j < 3; j++){
      b[ptId + j] = b[ptId + j] + (Ep[Id] * Bij[0][j]) + (Ep[Id + 1] * Bij[1][j]);
   }
}

/*-------------------------------------------------------------------------*/
