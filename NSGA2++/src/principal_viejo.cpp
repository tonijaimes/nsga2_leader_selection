/*
 * principal.cpp
 *
 *  Created on: 28/10/2009
 *      Author: antonio
 */

#include <cstdlib>
#include "RandUtils.h"
#include "NSGA.h"
#include "theMOPS.h"

typedef struct {
   int popsize;        // Population size.
   int ngen;           // Number of generations.
   int nreal;          // Number of variables using a real encoding.
   double pcross_real; // Crossover probability for real variables.
   double pmut_real;   // Mutation probability for real variables.
   double eta_c;
   double eta_m;
   int nbin;           // Number of variables using binary encoding.
   vector<int> nbits;  // Bits for each binary variable.
   double pcross_bin;  // Crossover probability for binary variables.
   double pmut_bin;    // Mutation probability for binary variables.
} NSGAParams;

void exitMessage(double value, const char *elementName);
void readParameters(int argc, char **argv, MOP *mop, Randomizer &r, NSGAParams &params);

int main(int argc, char **argv) {
   MOP *mop = new ZDT1();

   cout << "\nCreando Randomizer.\n";
   Randomizer r;
   NSGAParams p;

   readParameters(argc, argv, mop, r, p);

   cout << "\nInput data successfully entered, now performing initialization.\n";
   cout << "\nCreando NSGA.\n";
   NSGA nsga(mop, &r, p.popsize, p.ngen,
                      p.nreal, p.pcross_real, p.pmut_real, p.eta_c, p.eta_m,
                      p.nbin, p.pcross_bin, p.pmut_bin, p.nbits);
   cout << "\nNSGA creado.\n";

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

void readParameters(int argc, char **argv, MOP *mop, Randomizer &r, NSGAParams &params)
{
   if (argc<2)
   {
      cout << "\n Usage ./nsga2r random_seed \n";
      exit(1);
   }

   double seed = (double) atof(argv[1]);
   if (seed <= 0.0 || seed >= 1.0)
   {
      cout << "\n Entered seed value is wrong, seed value must be in (0,1) \n";
      exit(1);
   }

   r.setSeed(seed);

   cout << "\n Enter the problem relevant and algorithm relevant parameters...";
   cout << "\n Enter the population size (a multiple of 4) : ";
   scanf("%d", &(params.popsize));
   if (params.popsize < 4 || (params.popsize % 4) != 0)
      exitMessage(params.popsize, "population size");

   cout << "\n Enter the number of generations: ";
   scanf("%d", &(params.ngen));
   if ( params.ngen < 1 )
      exitMessage(params.ngen, "number of generations");

   cout << "\n There are a total of " << mop->getNumVariables()
   << " variables for problem " << mop->getName() << ".";

   cout << "\n Enter the number of real variables: ";
   scanf("%d", &(params.nreal));
   if (params.nreal < 0)
      exitMessage(params.nreal, "number of real variables");

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

   cout << "\n Enter the number of binary variables: ";
   scanf("%d", &(params.nbin));
   if ( params.nbin < 0 )
      exitMessage(params.nbin, "number of binary variables");

   if (params.nbin != 0)
   {
      params.nbits.resize(params.nbin);
      for (int i=0; i < params.nbin; ++i)
      {
         cout << "\n Enter the number of bits for binary variable: " << i+1;
         scanf ("%d", &(params.nbits[i]));
         if (params.nbits[i] < 1)
            exitMessage(params.nbits[i], "number of bits for binary variable");
      }
      cout << "\nEnter the probability of crossover of binary variable (0.6-1.0): ";
      scanf ("%lf", &(params.pcross_bin));
      if (params.pcross_bin < 0.0 || params.pcross_bin > 1.0)
         exitMessage(params.pcross_bin, "probability of crossover of binary variables");

      cout << "\n Enter the probability of mutation of binary variables (1/nbits): ";
      scanf ("%lf", &(params.pmut_bin));
      if (params.pmut_bin < 0.0 || params.pmut_bin > 1.0)
         exitMessage(params.pmut_bin, "probability  of mutation of binary variables");
   }

   if (params.nreal==0 && params.nbin==0) {
      cout << "\nNumber of real as well as binary variables, both are zero, hence exiting.\n";
      exit(1);
   }
}
