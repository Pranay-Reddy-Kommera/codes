/////////////////////////////////////////////////////////////////////////////
///   File           :          BA_io.cpp
///   Description    :
/////////////////////////////////////////////////////////////////////////////

#include "BA_declaration.h"

/*-------------------------------------------------------------------------*/
/*                       Declare external symbols                          */
/*-------------------------------------------------------------------------*/
unsigned int num_cam, num_pt, num_proj;
dtype *P, *projections, *calcProjections, *Pnew, *calcProjectionsnew;//, *projectionsFile;
unsigned int *ptidx, *camidx;//, *ptidxFile, *camidxFile;
dtype *Ep, *Epnew;
vector<dtype> errorVector, cgTime, cgIter;
dtype *U, *V, *b, *M;
dtype *valW, mu;
//unsigned int *rowPtrW;
bool stop;
dtype* delta;
dtype Epvector, Epvectornew, rho, vv;

/*-------------------------------------------------------------------------*/
/*          Function to read data from the text file provided.             */
/*-------------------------------------------------------------------------*/
bool readData(char* fileName)
{
   cout<<"Reading data from the given file ..."<<endl;
   
   ifstream infile(fileName);
   
//Read the first line in the file and obtain number of cameras, points and projections
   if(!(infile>>num_cam>>num_pt>>num_proj)){
      cout<<"Unable to read first line in specified format"<<endl;
      return 0;
   }
   cout<<"Number of camers: "<<num_cam<<"; points: "<<num_pt<<"; projections: "<<num_proj<<endl;

//Declare dimensions of the variables
   declareDim();

//Read complete data from file
   if(!(read(infile))){
      cout<<"Error reading data into declared dimensions"<<endl;
      return 0;
   }

   infile.close();

   return 1;
}

/*-------------------------------------------------------------------------*/
/*                  Declare Dimensions for the parameters                  */
/*-------------------------------------------------------------------------*/
void declareDim()
{
   //P = [(15*num_cam)+(3*num_pt)]
   P = new dtype[(15*num_cam)+(3*num_pt)];
   Pnew = new dtype[(15*num_cam)+(3*num_pt)];
   M = new dtype[(15*num_cam)+(3*num_pt)];

   //projections & calcProjections = [2*num_proj]
   projections = new dtype[2*num_proj];
//   projectionsFile = new dtype[2*num_proj];
   calcProjections = new dtype[2*num_proj];
   calcProjectionsnew = new dtype[2*num_proj];

   //ptidx = [num_proj]
   ptidx = new unsigned int[num_proj];
//   ptidxFile = new unsigned int[num_proj];

   //camidx = [num_proj]
   camidx = new unsigned int[num_proj];
//   camidxFile = new unsigned int[num_proj];

   //Ep vector = [2*num_proj]
   Ep = new dtype[2*num_proj];
   Epnew = new dtype[2*num_proj];

   //U vector = [num_cam * 225]
   U = new dtype[225 * num_cam];

   //V vector = [num_pt * 9]
   V = new dtype[9 * num_pt];

   //valW = [num_proj * 45]
   valW = new dtype[45 * num_proj];

   //rowPtrW = [num_cam + 1]
//   rowPtrW = new unsigned int[num_cam+1];

   //b = [(num_cam*15)+(num_pt*3)]
   b = new dtype[(num_cam*15)+(num_pt*3)];

   //delta = [(num_cam*15)+(num_pt*3)]
   delta = new dtype[(num_cam*15)+(num_pt*3)];
}

/*-------------------------------------------------------------------------*/
/*                         Read data into parameters                       */
/*-------------------------------------------------------------------------*/
bool read(ifstream& infile)
{
   if (!infile.is_open()) return 0;
   
   int idx = 0, projId = 0;;
   dtype x, y;
   int cam, pt;

   cout << "P vector is used to store camera parameters and 3D points." << endl;
   //Retrieve camidx, ptidx and projections
   for (int i = 0; i < num_proj; i++)
   {
      if (!(infile >> cam >> pt >> x >> y)){
         cout << "Could not read the camidx, ptidx and projections in the specified format" << endl;
         return 0;
      }

      projId = 2 * i;
      camidx[i] = cam; ptidx[i] = pt;
      projections[projId] = x; projections[projId + 1] = y;
   }

   //Retrieve camera parameters
   for (int i = 0; i < num_cam; i++)
   {
      dtype K[3], t[3], f, k1, k2;
      idx = (15 * i);

      if (!(infile >> K[0] >> K[1] >> K[2])){
         cout << "Could not read the Rodrigues vector in the specified format" << endl;
         return 0;
      }
      if (!(infile >> t[0] >> t[1] >> t[2])){
         cout << "Could not read the translation vector in the specified format" << endl;
         return 0;
      }
      if (!(infile >> f >> k1 >> k2)){
         cout << "Could not read the focal length and radial distortion parameters in the specified format" << endl;
         return 0;
      }

      setRotationMatrix(i, K, idx);
      P[idx + 9] = t[0]; P[idx + 10] = t[1]; P[idx + 11] = t[2];
      P[idx + 12] = f; P[idx + 13] = k1; P[idx + 14] = k2;
   }

   int id = 15 * num_cam;
   //Retrieve point coordinates
   for (int i = 0; i < num_pt; i++){
      dtype x, y, z;
      idx = id + (3 * i);
      infile >> x >> y >> z;

      P[idx] = x; P[idx + 1] = y; P[idx + 2] = z;
   }

   return 1;
}

/*-------------------------------------------------------------------------*/
/*                           Deallocate memory                             */
//-------------------------------------------------------------------------//
void clearData()
{
   //P=(a1', a2', a3',...,anum_cam',b1',b2',...,bnum_pt')
   //ai'=(r1, r2, r3, r4, r5, r6, r7, r8, r9, t1, t2, t3, focal, k1, k2)
   //bi'=(x, y, z)
   //P = [(15*num_cam)+(3*num_pt)]
   delete[] P;
   delete[] Pnew;
   delete[] M;

   //projections & calcProjections = [num_proj*2]
   delete[] projections;
//   delete[] projectionsFile
   delete[] calcProjections;
   delete[] calcProjectionsnew;

   //ptidx = [num_proj]
   delete[] ptidx;
//   delete[] ptidxFile;

   //camidx = [num_proj]
   delete[] camidx;
//   delete[] camidxFile;

   //Ep vector = [2*num_proj]
   delete[] Ep;
   delete[] Epnew;

   //U vector = [225 * num_cam]
   delete[] U;

   //V vector = [9 * num_pt]
   delete[] V;

   //valW = [num_proj * 45]
   delete[] valW;

   //rowPtrW = [num_cam + 1]
//   delete[] rowPtrW;

   //b = [(num_cam*15)+(num_pt*3)]
   delete[] b;

   //delta = [(num_cam*15)+(num_pt*3)]
   delete[] delta;
}

/*-------------------------------------------------------------------------*/
/*                     Write data to output text files                     */
//-------------------------------------------------------------------------//
bool writeInputData()
{
   char *end, *smallEnd; int projId=0;
   smallEnd = "----------";
   end = "-----------------------------------------";

   cout << "P vector is written into text documents." << endl;
   //Writing camera parameters to text file
   ofstream outfile("P_CameraParameters.txt");
   int idx;
   for (int i = 0; i < num_cam; i++){
      idx = 15 * i;
      outfile << end << "\n";
      outfile << smallEnd << "Camera " << i << smallEnd << "\n";
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
      outfile1 << smallEnd << "Point " << i << smallEnd << "\n";
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

   return 1;
}
