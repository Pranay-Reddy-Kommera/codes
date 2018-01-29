///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_dataRW.cpp
///   Description    :         Read and Write data.
//////////////////////////////////////////////////////////////////////////////////////

#include "BA_data.h"

/*********************************************************************/

void DeclareDim();
bool Read(ifstream& infile);

/*********************************************************************/

int num_cam, num_pt, num_proj;
dtype **R, **T, *focal, **rd, *points, *projections, *P, *calcProjections, *Pnew, *calcProjectionsnew, *Ep, *Epnew;
int *ptidx, *camidx, *ncpoints;
dtype *U, *V, *W, *g, mu, *delta, *Vinv;
dtype rho, Epvector, Epvectornew, vv;
bool stop;
vector<dtype> errorVector, rhoVector, muVector;

/*********************************************************************/
/*       Function to read data from the text file provided.          */
/*********************************************************************/
bool ReadData(char *fileName)
{
   //Read the text file and update parameters
   cout << "\nReading data from given file..." << endl;
   ifstream infile(fileName);
   //First read the file and obtain values for number of cameras, points and projections
   if (!(infile >> num_cam >> num_pt >> num_proj))
   {
      cout << "Unable to read the text file" << endl;
      return 0;
   }
   cout << "Number of cameras = " << num_cam << ", points = " << num_pt << " and projections = " << num_proj << endl;

   //First declare dimensions in all the arrays
   DeclareDim();

   //Read data from file
   if (!(Read(infile)))
   {
      cout << "Error reading data into declared dimensions" << endl;
      return 0;
   }

   infile.close();

   return 1;
}

/*********************************************************************/
/*  Write camera parameters, points and projections into text files. */
/*********************************************************************/
bool WriteInputData()
{
   char *end, *smallend; int projId = 0;
   smallend = "---------";
   end = "----------------------------------------------------";
   //cout << "\n";

#if defined _RTVector
   cout << "R, T, focal, rd are written into text documents." << endl;
   //Writing camera parameters to text file
   ofstream outfile("CameraParameters.txt");

   for (int i = 0; i < num_cam; i++){
      outfile << end << "\n";
      outfile << smallend<<"Camera " << i<<smallend<<"\n";
      outfile << "Rotation Matrix" << "\n";
      outfile << R[i][0] <<"\t"<< R[i][1] <<"\t"<< R[i][2]<<"\n";
      outfile << R[i][3] <<"\t"<< R[i][4] <<"\t"<< R[i][5]<<"\n";
      outfile << R[i][6] <<"\t"<< R[i][7] <<"\t"<< R[i][8]<<"\n";
      outfile << "Translation Vector" << "\n";
      outfile << T[i][0] << "\t" << T[i][1] << "\t" << T[i][2] << "\n";
      outfile << "focal = " <<focal[i]<< "\n";
      outfile << "Radial Distortion" << "\n";
      outfile << rd[i][0] << "\t" << rd[i][1]<<"\n";
   }

   outfile.close();

   //Writing point parameters to text file.
   ofstream outfile1("3D_Points.txt");

   for (int i = 0; i < num_pt; i++){
      outfile1 << smallend << "Point " << i << smallend << "\n";
      outfile1 << "x = " << points[i][0] << ", y = " << points[i][1] << ", z = " << points[i][2] << "\n";
   }

   outfile1.close();

   //Writing projection data to text file.
   ofstream outfile2("Projections.txt");

   for (int i = 0; i < num_proj; i++){
      projId=2*i;
      outfile2 << i << ". " << "camidx = " << camidx[i] << ", ptidx = " << ptidx[i] << ", projections: x = " << projections[projId] << ", y = " << projections[projId+1] << ".\n";
   }

   outfile2.close();
#endif

#if defined _PVector
   cout << "P vector is written into text documents." << endl;
   //Writing camera parameters to text file
   ofstream outfile("P_CameraParameters.txt");
   int idx;
   for (int i = 0; i < num_cam; i++){
      idx = 15 * i;
      outfile << end << "\n";
      outfile << smallend << "Camera " << i << smallend << "\n";
      outfile << "Rotation Matrix" << "\n";
      outfile << P[idx+0] << "\t" << P[idx+1] << "\t" << P[idx+2] << "\n";
      outfile << P[idx+3] << "\t" << P[idx+4] << "\t" << P[idx+5] << "\n";
      outfile << P[idx+6] << "\t" << P[idx+7] << "\t" << P[idx+8] << "\n";
      outfile << "Translation Vector" << "\n";
      outfile << P[idx+9] << "\t" << P[idx+10] << "\t" << P[idx+11] << "\n";
      outfile << "focal = " << P[idx+12] << "\n";
      outfile << "Radial Distortion" << "\n";
      outfile << P[idx+13] << "\t" << P[idx+14] << "\n";
   }

   outfile.close();

   ofstream outfile1("P_3DPoints.txt");

   int id = 15 * num_cam;
   for (int i = 0; i < num_pt; i++){
      idx = id + (3 * i);
      outfile1 << smallend << "Point " << i << smallend << "\n";
      outfile1 << "x = " << P[idx] << ", y = " << P[idx+1] << ", z = " << P[idx+2] << "\n";
   }

   outfile1.close();

   //Writing projection data to text file.
   ofstream outfile2("P_Projections.txt");

   for (int i = 0; i < num_proj; i++){
      projId = 2 * i;
      outfile2 << i << ". " << "camidx = " << camidx[i] << ", ptidx = " << ptidx[i] << ", projections: x = " << projections[projId] << ", y = " << projections[projId+1] << ".\n";
   }

   outfile2.close();
#endif

   return 1;
}

/*********************************************************************/
/*                    Read data into parameters                      */
/*********************************************************************/
bool Read(ifstream& infile)
{
   if (!infile.is_open()) return 0;
   int idx = 0, projId = 0;;
#if defined _RTVector
   cout << "R, T, focal and rd variables are used to store camera parameters." << endl;
   //Retrieve camidx, ptidx and projections
   for (int i = 0; i < num_proj; i++)
   {
      dtype x, y;
      int cam, pt;

      projId=2*i;
      infile >> cam >> pt >> x >> y;

      camidx[i] = cam; ptidx[i] = pt;
      projections[projId] = x; projections[projId+1] = y;
   }

   //Retrieve camera parameters
   for (int i = 0; i < num_cam; i++)
   {
      dtype K[3], t[3], f, k1, k2;
      infile >> K[0]; infile >> K[1]; infile >> K[2];
      infile >> t[0]; infile >> t[1]; infile >> t[2];
      infile >> f; infile >> k1; infile >> k2;

      SetRotationMatrix(i, K, idx);
      T[i][0] = t[0]; T[i][1] = t[1]; T[i][2] = t[2];
      focal[i] = f; rd[i][0] = k1; rd[i][1] = k2;
   }

   //Retrieve point coordinates
   for (int i = 0; i < num_pt; i++){
      dtype x, y, z;
      infile >> x >> y >> z;

      points[i][0] = x; points[i][1] = y; points[i][2] = z;
   }
#endif

#if defined _PVector
   cout << "P vector is used to store camera parameters and 3D points." << endl;
   //Retrieve camidx, ptidx and projections
   for (int i = 0; i < num_proj; i++)
   {
      dtype x, y;
      int cam, pt;

      infile >> cam >> pt >> x >> y;
      projId = 2 * i;
      camidx[i] = cam; ptidx[i] = pt;
      projections[projId] = x; projections[projId+1] = y;
   }

   //Retrieve camera parameters
   for (int i = 0; i < num_cam; i++)
   {
      dtype K[3], t[3], f, k1, k2;
      idx = (15 * i);
      infile >> K[0]; infile >> K[1]; infile >> K[2];
      infile >> t[0]; infile >> t[1]; infile >> t[2];
      infile >> f; infile >> k1; infile >> k2;

      SetRotationMatrix(i, K, idx);
      P[idx+9] = t[0]; P[idx+10] = t[1]; P[idx+11] = t[2];
      P[idx+12] = f; P[idx+13] = k1; P[idx+14] = k2;
   }

   int id = 15 * num_cam;
   //Retrieve point coordinates
   for (int i = 0; i < num_pt; i++){
      dtype x, y, z; 
      idx = id + (3 * i);
      infile >> x >> y >> z;

      P[idx] = x; P[idx+1] = y; P[idx+2] = z;
   }
#endif

   return 1;
}

/*********************************************************************/
/*                 Declare dimensions for parameters                 */
/*********************************************************************/
void DeclareDim()
{
#if defined _RTVector
   //R = [num_cam][9]
   //Memory for row elements
   R = new dtype *[num_cam];
   //Memory for column elements
   for (int i = 0; i < num_cam; i++){
      R[i] = new dtype[9];
   }

   //T = [num_cam][3]
   //Memory for row elements
   T = new dtype *[num_cam];
   //Memory for column elements
   for (int i = 0; i < num_cam; i++){
      T[i] = new dtype[3];
   }

   //focal = [num_cam]
   focal = new dtype[num_cam];

   //rd = [num_cam][2]
   //Memory for row elements
   rd = new dtype *[num_cam];
   //Memory for column elements
   for (int i = 0; i < num_cam; i++){
      rd[i] = new dtype[2];
   }

   //points = [num_pt][3]
   //Memory for row elements
   points = new dtype *[num_pt];
   //Memory for column elements
   for (int i = 0; i < num_pt; i++){
      points[i] = new dtype[3];
   }
#endif

#if defined _PVector
   //P=(a1', a2', a3',...,anum_cam',b1',b2',...,bnum_pt')
   //ai'=(r1, r2, r3, r4, r5, r6, r7, r8, r9, t1, t2, t3, focal, k1, k2)
   //bi'=(x, y, z)
   //P = [(15*num_cam)+(3*num_pt)]
   P = new dtype[(15 * num_cam) + (3 * num_pt)];
   Pnew = new dtype[(15 * num_cam) + (3 * num_pt)];
#endif

   //projections & calcProjections = [num_proj*2]
   projections = new dtype [2*num_proj];
   calcProjections = new dtype[2 * num_proj];
   calcProjectionsnew = new dtype[2 * num_proj];
   Ep = new dtype[2 * num_proj];
   Epnew = new dtype[2 * num_proj];

   //ptidx = [num_proj]
   ptidx = new int [num_proj];

   //camidx = [num_proj]
   camidx = new int [num_proj];

   //ncpoints = [num_pt]
   ncpoints = new int[num_pt];

   //U vector
   U = new dtype[num_cam*225];

   //V vector
   V = new dtype[num_pt*9];
   Vinv = new dtype[num_pt*9];

   //W vector
   W = new dtype[num_proj*45];

   //g vector
   g = new dtype[(num_cam*15)+(num_pt*3)];

   //delta vector
   delta = new dtype[(num_cam * 15) + (num_pt * 3)];

}

void ClearData()
{
#if defined _RTVector
   //R = [num_cam][9]
   for (int i = 0; i < num_cam; i++){
      delete[] R[i];
   }
   delete[] R;

   //T = [num_cam][3]
   for (int i = 0; i < num_cam; i++){
      delete[] T[i];
   }
   delete[] T;

   //focal = [num_cam]
   delete[] focal;

   //rd = [num_cam][2]
   for (int i = 0; i < num_cam; i++){
      delete[] rd[i];
   }
   delete[] rd;

   //points = [num_pt][3]
   for (int i = 0; i < num_pt; i++){
      delete[] points[i];
   }
   delete[] points;
#endif

#if defined _PVector
   //P=(a1', a2', a3',...,anum_cam',b1',b2',...,bnum_pt')
   //ai'=(r1, r2, r3, r4, r5, r6, r7, r8, r9, t1, t2, t3, focal, k1, k2)
   //bi'=(x, y, z)
   //P = [(15*num_cam)+(3*num_pt)]
   delete[] P;
   delete[] Pnew;
#endif

   //projections & calcProjections = [num_proj*2]
   delete[] projections;
   delete[] calcProjections;
   delete[] calcProjectionsnew;
   delete[] Ep;
   delete[] Epnew;

   //ptidx = [num_proj]
   delete[] ptidx;

   //camidx = [num_proj]
   delete[] camidx;

   //ncpoints = [num_pt]
   delete[] ncpoints;

   //U vector
   delete[] U;

   //V vector
   delete[] V;
   delete[] Vinv;

   //W vector
   delete[] W;

   //g vector
   delete[] g;

   //delta vector
   delete[] delta;
}
