/*
 * principal.cpp
 *
 *  Created on: 28/10/2009
 *      Author: antonio
 */

#include <NSGA.h>
#include <cstdlib>
#include "RandUtils.h"
#include "SensorNetwork.h"

//#include "MOPsCatalog.h"
//#include "theMOPS.h"

//typedef struct {
//   int popsize;        // Population size.
//   int ngen;           // Number of generations.
//   int nreal;          // Number of variables using a real encoding.
//   double pcross_real; // Crossover probability for real variables.
//   double pmut_real;   // Mutation probability for real variables.
//   double eta_c;
//   double eta_m;
//   int nbin;           // Number of variables using binary encoding.
//   vector<int> nbits;  // Bits for each binary variable.
//   double pcross_bin;  // Crossover probability for binary variables.
//   double pmut_bin;    // Mutation probability for binary variables.
//} NSGAParams;

void exitMessage(double value, const char *elementName);
void readRandomSeed(int argc, char **argv, double &seed);
void readMOPParameters(int argc, char **argv, MOP **mop);
void readNSGAParameters(int argc, char **argv, MOP *mop, NSGAParams &params);

int main(int argc, char **argv) {
   MOP *mop = new SensorNet(200, 2000, 2000);

   cerr << "\nReading parameters..." << endl;
   if (argc < 2) {
      cout << "\n Usage ./nsga2r random_seed\n";
      exit(1);
   }

   Randomizer r;
   double seed;
   readRandomSeed(argc, argv, seed);
   r.setSeed(seed);

   readMOPParameters(argc, argv, &mop);

   NSGAParams p;
   readNSGAParameters(argc, argv, mop, p);

   cout << "\nInput data successfully entered, now performing initialization.\n";
   NSGA nsga(mop, &r, p);

   cout << "\nNSGA-II was successfully created.\n";

   cout << "\nNow performing optimization...\n";
   nsga.optimize();

   cout << "\nGenerations finished, now reporting solutions.";
   cout << "\nRoutine successfully exited \n";

   return 0;
}

void exitMessage(double value, const char *elementName)
{
   cerr << "\n" << elementName << " read is: " << value;
   cout << "\nWrong value for " << elementName << ", hence exiting.\n";
   exit (1);
}

void readRandomSeed(int argc, char **argv, double &seed) {
   seed = (double) atof(argv[1]);
   if (seed <= 0.0 || seed >= 1.0) {
      cout << "\n Entered seed value is wrong, seed value must be in (0,1) \n";
      exit(0);
   }
}


void readMOPParameters(int argc, char **argv, MOP **mop)
{
   cout << "\nThe MOP selected is " << (*mop)->getName() << endl;

//   *mop = mops.selectMOP(mopName);
//   if (*mop == NULL) {
//      cout << "MOP '" << mopName << "' does not exists. Available MOPs: \n\n";
//      mops.displayMOPList();
//      exit(0);
//   }
//
//   //TODO: check if there are more parameters for MOP, namely # variables and #objectives.
//   if ((*mop)->isScalableObjSpace()) {
//      if (argc >= 6) { //Read num. of variables and objectives.
//         int nObjectives = (int) atoi(argv[4]);
//         int nVariables  = (int) atoi(argv[5]);
//
//         (*mop)->setNumVariables(nVariables);
//         (*mop)->setNumObjectives(nObjectives);
//      }
//   }
}

void readNSGAParameters(int argc, char **argv, MOP *mop, NSGAParams &params)
{
   cout << "\nEnter the problem relevant and algorithm relevant parameters...";
   cout << "\nEnter the population size (a multiple of 4) : ";
   scanf("%d", &(params.popsize));
   if (params.popsize < 4 || (params.popsize % 4) != 0)
      exitMessage(params.popsize, "population size");

   cout << "\n Enter the number of generations: ";
   scanf("%d", &(params.ngen));
   if ( params.ngen < 1 )
      exitMessage(params.ngen, "number of generations");

   cout << "\n There is a total of " << mop->getNumVariables()
   << " VARIABLES for problem " << mop->getName() << ".";

   cout << "\n Enter the number of real variables: ";
   //scanf("%d", &(params.nreal));

   params.nreal = 0;
   if (params.nreal < 0 || params.nreal > mop->getNumVariables())
      exitMessage(params.nreal, "number of real variables");

   params.nbin = mop->getNumVariables() - params.nreal;
   if ( params.nbin < 0 ) // Esto no debe ser necesario.
      exitMessage(params.nbin, "number of binary variables");

   cout << "\n The " << params.nbin << " remaining variables will have binary representation.\n";

   if (params.nreal != 0)
   {
      cout << "\n Enter the probability of crossover of real variable (0.6-1.0): ";
      scanf ("%lf", &(params.pcross_real));
      if (params.pcross_real < 0.0 || params.pcross_real > 1.0)
         exitMessage(params.pcross_real, "probability of crossover of real variables");

      cout << "\n Enter the probability of mutation of real variables (1/nreal): ";
      scanf ("%lf", &(params.pmut_real));
      if (params.pmut_real < 0.0 || params.pmut_real > 1.0)
         exitMessage(params.pmut_real, "probability of mutation of real variables");

      cout << "\n Enter the value of distribution index for crossover (5-20): ";
      scanf ("%lf", &(params.eta_c));
      if (params.eta_c <= 0)
         exitMessage(params.eta_c, "distribution index for crossover");

      cout << "\n Enter the value of distribution index for mutation (5-50): ";
      scanf ("%lf", &(params.eta_m));
      if (params.eta_m <= 0)
         exitMessage(params.eta_m, "distribution index for mutation");
   }

   if (params.nbin != 0)
   {
      //ESTO ES PARA UNA VARIABLE BINARIA CON UN BIT POR CADA SENSOR
      params.nbits.resize(1, params.nbin);   // vector nbits de longitud 1 con el valor de nbin (Núm sensores).
      params.nbin = 1; // Suponer que es una sola variable
/*
      params.nbits.resize(params.nbin);
      for (int i=0; i < params.nbin; ++i)
      {
         cout << "\n Enter the number of bits for binary variable " << i+1 << ": ";
         scanf ("%d", &(params.nbits[i]));
         if (params.nbits[i] < 1)
            exitMessage(params.nbits[i], "number of bits for binary variable");
      }
      */
      cout << "\nEnter the probability of crossover of binary variable (0.6-1.0): ";
      scanf ("%lf", &(params.pcross_bin));
      if (params.pcross_bin < 0.0 || params.pcross_bin > 1.0)
         exitMessage(params.pcross_bin, "probability of crossover of binary variables");

      cout << "\n Enter the probability of mutation of binary variables (1/nbits): ";
      params.pmut_bin = 1.0 / params.nbits[0];
      //scanf ("%lf", &(params.pmut_bin));
      if (params.pmut_bin < 0.0 || params.pmut_bin > 1.0)
         exitMessage(params.pmut_bin, "probability  of mutation of binary variables");
   }
}
