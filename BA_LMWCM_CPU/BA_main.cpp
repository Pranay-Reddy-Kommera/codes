///////////////////////////////////////////////////////////////////////////////////////
///   File           :         BA_main.cpp
///   Description    :         Bundle Adjustment implementation w.r.t. SBA
//////////////////////////////////////////////////////////////////////////////////////

#include "BA_data.h"

/*********************************************************************/

int main(int argc, char* argv[])
{
   //auto start = chrono::high_resolution_clock::now();
   cout << "Bundle Adjustment has started ... " << endl;
   int itr = 0;

   //Read the file name from the command line
   char* fileName;
   fileName = argv[1];
   //fileName = "BA_257.txt"; //For now we hardcode the file name
   //fileName = "BA_21.txt"; //For now we hardcode the file name

   //Read data from the file
   if (!(ReadData(fileName)))
   {
      cout << "Error retrieving data from " << fileName << endl;
      //goto CodeEnd;
      return 0;
   }

   //Write data into text file for confirmation
#if defined _Write
   cout << "\nInput is written into text documents." << endl;
   if (!(WriteInputData()))
   {
      cout << "Error writing data to " << fileName << endl;
      //goto CodeEnd;
      return 0;
   }
#endif

   cout<<"Data is read into memory"<<endl;

   //Calculate measured projections
   //Calculating f(P)
   MeasuredProjections();

   //Calculate RMS Error
   //Calculate Ep
   errorCalculation();
   RMSErrorCalculation();

   //Calculate JtJ and JtEp
   UVWgCalculation();

   //Stopping criteria
   stop = stopCriteria();
   //cout << "stop is " << stop << endl;

   //Calculate mu
   mu = muCalculation();
   //cout << "mu is " << mu << endl;

   //Calculate number of camers each points are in
   numCamPtsCalculations();
   cout<<"Starting iterations..."<<endl;

//************   Iteration starts here   **************
   auto start = chrono::high_resolution_clock::now();

   while (!stop && itr < kmax){
      itr = itr + 1;
      while (1){
         //cout<<"Inner while loop"<<endl;
         //Calculate delta=[deltaA' deltaB']

         deltaCalculation();

         //cout<<"Delta is calculated"<<endl;

         if (vecNorm(delta, ((num_cam * 15) + (num_pt * 3))) <= (e3*(vecNorm(P, ((num_cam * 15) + (num_pt * 3))) + e3))){
            stop = 1;
         }
         else{
            //Update Parameters Pnew with P
            updateParametersPnew();

            //Calculate measured projections
            MeasuredProjectionsnew();

            //Calculate Epnew
            EpnewCalculation();

            //rho calculation
            Epvector = vecNorm(Ep, (2 * num_proj)); Epvectornew = vecNorm(Epnew, (2 * num_proj));
            rho = ((Epvector*Epvector) - (Epvectornew*Epvectornew)) / (calculateDenom());
            //cout << "rho = " << rho <<"and denom is "<<calculateDenom()<< endl;
            //cout << "Epvector is " << Epvector << "and Epvectornew is " << Epvectornew << endl;
            //rhoVector.push_back(rho);

            if (rho > 0){
               //cout<<"Entered as rho > 0"<<endl;
               stop = ((Epvector - Epvectornew) < (e2*Epvector));
               //stop = ((vecNorm(Ep, (2 * num_proj)) - vecNorm(Epnew, (2 * num_proj))) < (e2*vecNorm(Ep, (2 * num_proj))));
					
               //Update Parameters P with Pnew
               updateParametersP();

               ////Update CalProjections with CalProjectionsnew
               //updateCalcProjections();

               //Calculate measured projections
               MeasuredProjections();

               ////Update Ep with Epnew
               //updateEp();

               //Calculate RMS error
               errorCalculation();
               RMSErrorCalculation();

               //Calculate JtJ and JtEp
               UVWgCalculation();

               //Determine stopping criteria
               stop = (stop) || (stopCriteria());
               if ((1.0 / 3) >= (1 - ((2 * rho - 1)*(2 * rho - 1)*(2 * rho - 1)))){
                  mu = mu*(1.0 / 3); vv = 2;
               }
               else{
                  mu = mu*(1 - ((2 * rho - 1)*(2 * rho - 1)*(2 * rho - 1))); vv = 2;
               }
            }
            else{
               mu = mu*vv; vv = 2 * vv;
            }
         }
         if (rho > 0 || stop) {
            break;
         }
      }
      stop = (vecNorm(Ep, (2 * num_proj)) <= e1);
   }

   auto end = chrono::high_resolution_clock::now();
   cout << "Time taken for "<<itr<< " iterations in milliseconds is " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;

//************   Iteration ends here   ******************
   //Print errorvector
   cout << "The error vector is " << endl;
   for (int l = 0; l < errorVector.size(); l++){
      cout << errorVector[l] << endl;
   }

   //Clear the data
   ClearData();
	
   cout << "\nBundle Adjustment has been completed ... " << endl;
   //auto end = chrono::high_resolution_clock::now();
   //auto elapsed = chrono::duration_cast < chrono::milliseconds > (end - start);
   //cout << "Time taken in milliseconds is " <<elapsed.count()<< endl;
   //cout << "Time taken in milliseconds is " << chrono::duration_cast<chrono::milliseconds>(end-start).count() << endl;
}

