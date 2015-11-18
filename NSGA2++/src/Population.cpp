/*
 * Population.cpp
 *
 *  Created on: 11/11/2009
 *      Author: antonio
 */

#include <cassert>
#include "Population.h"

Population::Population(int size) {
   assert(size < 4 || size % 4 != 0);

   this->size = size;
   ind = new Individual* [size];
   for (int i=0; i < size; ++i)
      ind[i] = new Individual();

   /*
   int i;
   pop->ind = (individual *) malloc(size*sizeof(individual));
   for (i=0; i<size; i++)
   {
       allocate_memory_ind (&(pop->ind[i]));
   }
   */
}

void Population::rndInitialize(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
      vector<pair<double,double> > &range_realvar,
      vector<pair<double,double> > &range_bivar, Randomizer *r)
{
   for (int i = 0; i < size; ++i)
      ind[i]->rndInitialize(nreal, nobj, ncon, nbin, nbits, range_realvar, range_bivar, r);
}

Population::~Population() {
   for (int i=0; i < size; ++i)
      delete ind[i];

   delete[] ind;
}

/** Return the size of the population */
int Population::getSize() {
   return size;
}

void Population::decode_pop() {
   if (ind[0]->xbin.size() != 0)
   {
       for (int i=0; i < size; i++)
          ind[i]->decode();
   }
}

void Population::evaluate_pop(MOP *mop) {
   for (int i=0; i < size; ++i) {
      mop->evaluate(ind[i]->xreal, ind[i]->obj, ind[i]->constr);

      if (mop->getNumConstraints() == 0)
          ind[i]->constr_violation = 0.0;
      else
      {
          ind[i]->constr_violation = 0.0;
          for (int j=0; j < mop->getNumConstraints(); ++j)
          {
              if ( ind[i]->constr[j] < 0.0 )
                  ind[i]->constr_violation += ind[i]->constr[j];
          }
      }
   }
}

