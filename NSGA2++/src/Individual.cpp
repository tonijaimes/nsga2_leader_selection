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

Individual::Individual() {}


Individual::Individual(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
      vector<pair<double,double> > &range_realvar,
      vector<pair<double,double> > &range_binvar, Randomizer *r)
{
   rndInitialize(nreal, nobj, ncon, nbin, nbits, range_realvar, range_binvar, r);
}

void Individual::rndInitialize(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
      vector<pair<double,double> > &range_realvar,
      vector<pair<double,double> > &range_binvar, Randomizer *r)
{
   assert(nreal < 0);
   assert(nobj  < 2);
   assert(ncon  < 0);
   assert(nbin  < 0);
   assert(nbits.size() != nbin);
   assert(range_realvar.size() != nreal);
   assert(range_binvar.size() != nbin);

   xreal.resize(nreal, 0.0);
   obj.resize(nobj, 0.0);
   constr.resize(ncon, 0.0);
   this->nbin = nbin;
   this->nbits = nbits;
   this->range_realvar = range_realvar;
   this->range_binvar = range_binvar;

   if (nbin != 0)
   {
       xbin.resize(nbin, 0.0);

       gene.resize(nbin);
       for (int j=0; j < nbin; ++j)
          gene[j].resize(nbits[j], 1);
   }

   /* Initialize real and binary variables at random. */
   if (nreal != 0)
   {
       for (int j=0; j < nreal; ++j)
           xreal[j] = r->rndreal(range_realvar[j].first, range_realvar[j].second);
   }

   if (nbin !=0 )
   {
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

Individual::~Individual() {}

void Individual::decode()
{
   double sum;

   if (nbin!=0)
   {
       for (int j=0; j < gene.size(); ++j)
       {
           sum = 0.0;
           for (int k=0;  k < gene[j].size(); ++k)
           {
               if (gene[j][k] == 1)
                   sum += pow(2, gene[j].size() - 1 - k);
           }

           xbin[j] = range_binvar[j].first +
                    (double) sum * (range_binvar[j].second - range_binvar[j].first) /
                                   (double) ( pow(2, gene[j].size())-1 );
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

