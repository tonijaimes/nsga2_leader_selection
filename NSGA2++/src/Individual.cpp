/*
 * Individual.cpp
 *
 *  Created on: 11/11/2009
 *      Author: antonio
 */

#include "Individual.h"
#include "RandUtils.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

Individual::Individual() {}


Individual::Individual(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
      vector<pair<double,double> > &range_realvar,
      vector<pair<double,double> > &range_binvar, Randomizer *r)
{
   rndInitialize(nreal, nobj, ncon, nbin, nbits, range_realvar, range_binvar, r);
}

void Individual::cloneInitialize(Individual &ind) {
   assert(ind.xreal.size() >= 0);
   assert(ind.xbin.size()  >= 0);
   assert((ind.xreal.size() + ind.xbin.size())  >= 1); /* at least one variable. */
   assert(ind.obj.size()  >= 2);
   assert(ind.constr.size()  >= 0);
   assert(ind.nbits.size() == (unsigned) ind.xbin.size());
   assert(ind.range_realvar.size() == (unsigned) ind.xreal.size());
//   assert(ind.range_binvar.size() == (unsigned) ind.xbin.size());

   this->nbin = ind.xbin.size();
   this->nbits = ind.nbits;
   this->range_realvar = ind.range_realvar;
   this->range_binvar = ind.range_binvar;

   xreal  = ind.xreal;
   xbin   = ind.xbin;
   obj    = ind.obj;
   constr = ind.constr;
   gene   = ind.gene;
}

void Individual::rndInitialize(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
      vector<pair<double,double> > &range_realvar,
      vector<pair<double,double> > &range_binvar, Randomizer *r)
{
//	cout << "range_binvar.size()=" << range_binvar.size() << "\n"
//	     << "nbin= " << nbin << endl;

   assert(nreal >= 0);
   assert(nbin  >= 0);
   assert((nreal + nbin)  >= 1);
   assert(nobj  >= 2);
   assert(ncon  >= 0);
   assert(nbits.size() == (unsigned)nbin);
   assert(range_realvar.size() == (unsigned)nreal);
//   assert(range_binvar.size() == (unsigned)nbin);

   xreal.resize(nreal, 0.0);
   obj.resize(nobj, 0.0);
   constr.resize(ncon, 0.0);
   this->nbin = nbin;
   this->nbits = nbits;
   this->range_realvar = range_realvar;
   this->range_binvar = range_binvar;

   //cout << "Range real size= "<< this->range_realvar.size() << " nreal=" << nreal;
   //cout << "Range bin size= "<< this->range_binvar.size() << " nbin=" << this->nbin;

   if (nbin != 0) {
       xbin.resize(nbin, 0.0);
       gene.resize(nbin);
       for (int j=0; j < nbin; ++j)
          gene[j].resize(nbits[j], 1);
   }

   /* Initialize real and binary variables at random. */
   //cout << "Random seed: " << r->getSeed() << endl;
   double x;
   if (nreal != 0) {
       for (int j=0; j < nreal; ++j) {
           x = r->rndreal(range_realvar[j].first, range_realvar[j].second);
           xreal[j] = x;
           //xreal[j] = r->rndreal(range_realvar[j].first, range_realvar[j].second);
           //cout << x << "\t";
       }
       //cout << "\n";
   }

   if (nbin !=0 ) {
       for (int j=0; j < nbin; ++j)
           for (int k=0; k < nbits[j]; ++k)
               gene[j][k] = (r->randomperc() <= 0.5) ? 0 : 1;
   }

   /*   void NSGA::allocate_memory_ind (individual *ind) {
      int j;
      if (nreal != 0)
      {
          ind->xreal = (double *)malloc(nreal*sizeof(double));
      }
      if (nbin != 0)
      {
          ind->xbin = (double *)malloc(nbin*sizeof(double));
          ind->gene = (int **)malloc(nbin*sizeof(int));
          for (j=0; j<nbin; j++)
          {
              ind->gene[j] = (int *)malloc(nbits[j]*sizeof(int));
          }
      }
      ind->obj = (double *)malloc(nobj*sizeof(double));
      if (ncon != 0)
      {
          ind->constr = (double *)malloc(ncon*sizeof(double));
      }
      return;
   }
*/
/*
   int j, k;
   if (nreal!=0)
   {
       for (j=0; j<nreal; j++)
       {
           ind->xreal[j] = r->rndreal (min_realvar[j], max_realvar[j]);
       }
   }
   if (nbin!=0)
   {
       for (j=0; j<nbin; j++)
       {
           for (k=0; k<nbits[j]; k++)
           {
               if (r->randomperc() <= 0.5)
               {
                   ind->gene[j][k] = 0;
               }
               else
               {
                   ind->gene[j][k] = 1;
               }
           }
       }
   }
   return;
*/
}


void Individual::initialize(vector<double>& x0,
                            int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
                            vector<pair<double,double> > &range_realvar,
                            vector<pair<double,double> > &range_binvar)
{
// cout << "range_binvar.size()=" << range_binvar.size() << "\n"
//      << "nbin= " << nbin << endl;

   //This is only for real enconded variables, so nbin must be 0.
   assert(nreal >= 1 && (int)x0.size() == nreal);
   assert(nbin  == 0);
   assert(nobj  >= 2);
   assert(ncon  >= 0);
   assert(nbits.size() == (unsigned)nbin);
   assert(range_realvar.size() == (unsigned)nreal);
//   assert(range_binvar.size() == (unsigned)nbin);

   xreal.resize(nreal, 0.0);
   obj.resize(nobj, 0.0);
   constr.resize(ncon, 0.0);
   this->nbin = nbin;
   this->nbits = nbits;
   this->range_realvar = range_realvar;
   this->range_binvar = range_binvar;

   if (nreal != 0) {
       for (int j = 0; j < nreal; ++j) {
           assert(x0[j] >= range_realvar[j].first && x0[j] <= range_realvar[j].second);
           xreal[j] = x0[j];
       }
   }
}

Individual::~Individual() {}

void Individual::printObjs(ostream &output) {
	for (unsigned i = 0; i < obj.size(); ++i)
		output << obj[i] << " ";
}

void Individual::copyIndividual(Individual &ind) {
   rank = ind.rank;
   constr_violation = ind.constr_violation;
   crowd_dist = ind.crowd_dist;
   obj = ind.obj;
   constr = ind.constr;
   xreal = ind.xreal;
   xbin = ind.xbin;
   nbin = ind.nbin;
   gene = ind.gene;
   nbits = ind.nbits;
   range_realvar = ind.range_realvar;
   range_binvar = ind.range_binvar;
}

void Individual::decode()
{
   double sum;

   if (nbin!=0)
   {
       for (unsigned j=0; j < gene.size(); ++j)
       {
           sum = 0.0;
           for (unsigned k=0;  k < gene[j].size(); ++k)
           {
               if (gene[j][k] == 1)
                   sum += pow(2, gene[j].size() - 1 - k);
           }

           xbin[j] = range_binvar[j].first +
                    (double) sum * (range_binvar[j].first, range_binvar[j].second) /
                                   (double) (pow(2, gene[j].size()) - 1);
       }
   }
}

void Individual::setRank(double rank) {
   this->rank = rank;
}

int Individual::getRank() {
   return rank;
}

void Individual::setCrowd_dist(double crowd_dist) {
   this->crowd_dist = crowd_dist;
}

double Individual::getCrowd_dist() {
   return crowd_dist;
}

void Individual::setConstr_violation(double constr_viol) {
   this->constr_violation = constr_viol;
}

double Individual::getConstr_violation() {
   return constr_violation;
}

