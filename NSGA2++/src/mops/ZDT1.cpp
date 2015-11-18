#include "ZDT1.h"
#include <cmath>
#include <iostream>

using namespace std;

// ZDT1 with a number of variables defined by the user (30 by default),
// and 2 objectives and 0 constraints.
ZDT1::ZDT1(int numVars) : MOP("ZDT1", numVars, 2, 0)
{
   // 0 <= xi <= 1, for i= 1...numVars
   xRanges.assign(nVariables, make_pair(0.0, 1.0));
}

ZDT1::~ZDT1() {}

void ZDT1::evaluate(double const *x, double *eval, double *gcons) const
{
   double f1, f2, g, h;

   f1 = x[0];

   g = 0.0;
   for (int i = 1; i < nVariables; ++i)
       g += x[i];

   g = 1.0 + 9.0 * (g/(nVariables-1));
   h = 1.0 - sqrt(f1/g);
   
   f2 = g * h;
   
   eval[0] = f1;
   eval[1] = f2;  
   
//   cout << "\nLlamando ZDT1 con arreglos." << endl;
}

void ZDT1::evaluate(vector<double> const &x, vector<double> &eval,  vector<double> &gcons) const
{
   double f1, f2, g, h;

   f1 = x[0];

   g = 0.0;
   for (int i = 1; i < nVariables; ++i)
       g += x[i];

   g = 1.0 + 9.0 * (g/(nVariables-1));
   h = 1.0 - sqrt(f1/g);
   
   f2 = g * h;
   
   eval[0] = f1;
   eval[1] = f2;  
   
  // cout << "\nLlamando ZDT1 con vectores." << endl;   
}

