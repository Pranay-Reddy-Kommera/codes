/////////////////////////////////////////////////////////////////////////////
///   File           :          BA_main.cpp
///   Description    :
/////////////////////////////////////////////////////////////////////////////

#include "BA_declaration.h"

int main(int argc, char* argv[])
{
   cout<<"-------------------------------------------"<<endl;
   cout<<"Bundle Adjustment started..."<<endl;
   cout<<"-------------------------------------------\n"<<endl;
   cout<<"Algorithm: Conjugate Gradient"<<endl;
   cout<<"Preconditioner: NO"<<endl;

   char *fileName=argv[1];
   cout<<"Input File: "<<fileName<<endl;
   cout<<"-------------------------------------------\n"<<endl;

/*-------------------------------------------------------------------------*/
/*                         Read data into variables                        */
/*-------------------------------------------------------------------------*/
   //Start time stamp
   auto start = chrono::high_resolution_clock::now();

   //Read data into variables from input file.
   if(!(readData(fileName))){
      cout<<"Error retrieving data from "<<fileName<<endl;
      return 0;
   }
   //Stop time
   auto end = chrono::high_resolution_clock::now();
   cout<<"Time taken to read data is "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" milliseconds"<<endl;
   cout<<"Data retrieving from the given file is successful"<<endl;
   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*        Convert camidx, ptidx and projections into cam format            */
/*-------------------------------------------------------------------------*/
//   changeFormat();

/*-------------------------------------------------------------------------*/
/*                        Write data into variables                        */
/*-------------------------------------------------------------------------*/
#if defined _write
   //Start time stamp
   start = chrono::high_resolution_clock::now(); 

   cout<<"Input is written into text documents"<<endl;
   if(!(writeInputData())){
      cout<<"Error writing data to text documents"<<endl;
      return 0;
   }

   //Stop time
   end = chrono::high_resolution_clock::now();
   cout<<"Time taken to write data is "<<chrono::duration_cast<chrono::milliseconds>(end-start).count()<<" milliseconds"<<endl;

   cout<<"-------------------------------------------"<<endl;
#endif

/*-------------------------------------------------------------------------*/
/*                    Calculate measured projections                       */
/*-------------------------------------------------------------------------*/
   //Calculate the calcProjections from camera and point parameters
   //Calculate f(P)
   measuredProjections();
   cout<<"\nMeasured projections are calculated"<<endl;
   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*               Calculate error vector and rms values                     */
/*-------------------------------------------------------------------------*/
   //Calculate error vector and rms value of error
   errorVectorCalculation();
   errorRMSCalculation();
   cout<<"\nError vector and rms value of error is calculated"<<endl;
   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*                     Calculate A matrix and b vector                     */
/*-------------------------------------------------------------------------*/
   //A and b
   calculateAb();
   cout<<"\nMatrix A and vector b are calculated"<<endl;
   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*                     Calculate stopping criteria                         */
/*-------------------------------------------------------------------------*/
   stop = stopCriteria();
   cout<<"\nStopping criteria is calculated: "<<stop<<endl;
   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*                     Calculating damping term mu                         */
/*-------------------------------------------------------------------------*/
   mu = muCalculation();
   cout<<"\nDamping term is calculated: "<<mu<<endl;
   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*                              Augment U and V                            */
/*-------------------------------------------------------------------------*/
//   augmentUV();
//   cout<<"\nU and V matrices are augmented"<<endl;
//   cout<<"-------------------------------------------"<<endl;

/*-------------------------------------------------------------------------*/
/*                     While loop to run the algorithm                     */
/*-------------------------------------------------------------------------*/
   int itr = 0;
   int size = (num_cam*15)+(num_pt*3);

   start = chrono::high_resolution_clock::now();

   while(!stop && itr < kmax){
      itr = itr + 1;
      while(1){

         augmentUV();
         //Compute jacobi preconditioner M^-1
         computeJacobiPreconditioner();
         //Conjugate gradient without preconditioner
         cgwjp();

//cout<<"vector delta "<<vecNorm(delta, size)<<endl;
//cout<<"vector P "<<(e1*(vecNorm(P,size)+e1))<<endl;

         if (vecNorm(delta, size) <= (e1*(vecNorm(P, size) + e1))){
            stop = 1;
//cout<<"*** Executed ***"<<endl;
         }
         else{
            //Update Parameters Pnew with P
            updateParametersPnew();

            //Calculate measured projections
            measuredProjectionsnew();

            //Calculate Epnew
            EpnewCalculation();

            //rho calculation
            Epvector = vecNorm(Ep, (2 * num_proj)); Epvectornew = vecNorm(Epnew, (2 * num_proj));
            rho = ((Epvector*Epvector) - (Epvectornew*Epvectornew)) / (calculateDenom());
///				cout << "rho = " << rho <<"and denom is "<<calculateDenom()<< endl;
///				cout << "Epvector is " << Epvector << "and Epvectornew is " << Epvectornew << endl;
///				rhoVector.push_back(rho);

				
            if (rho > 0){
               //cout << "Inside rho loop" << endl;
               stop = ((Epvector - Epvectornew) < (e2*Epvector));
///               stop = ((vecNorm(Ep, (2 * num_proj)) - vecNorm(Epnew, (2 * num_proj))) < (e2*vecNorm(Ep, (2 * num_proj))));
					
               //Update Parameters P with Pnew
               updateParametersP();

///            //Update CalProjections with CalProjectionsnew
///            updateCalcProjections();

               //Calculate measured projections
               measuredProjections();

///            //Update Ep with Epnew
///            updateEp();

               //Calculate RMS error
               errorVectorCalculation();
               errorRMSCalculation();

               //Calculate A and b
               calculateAb();

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

   end = chrono::high_resolution_clock::now();
   cout << "Time taken for "<<itr<< " iterations in milliseconds is " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << endl;
/*-------------------------------------------------------------------------*/

   printCGDetails();
   printErrorVector();
/*-------------------------------------------------------------------------*/
   cout<<"-------------------------------------------"<<endl;
   //Free memory
   clearData();
   cout<<"\nMemory deallocated"<<endl;
   cout<<"-------------------------------------------"<<endl;

   cout<<"-------------------------------------------"<<endl;
   cout<<"Bundle Adjustment ended..."<<endl;
   cout<<"-------------------------------------------"<<endl;

   return 1;
}

