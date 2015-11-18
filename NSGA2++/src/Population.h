/*
 * Poblacion.h
 *
 *  Created on: 28/10/2009
 *      Author: antonio
 */

#ifndef POBLACION_H_
#define POBLACION_H_

#include "theMOPS.h"
#include "Individual.h"
#include "RandUtils.h"

/*
typedef struct
{
    individual *ind;
    int size;
}
population;
*/

class Population {
public:
   int size;

   int nreal;
   int nobj;
   int ncon;
   int nbin;
   vector<int> nbits;
   vector<pair<double,double> >  range_realvar;

public:
   Individual **ind;

   Population(int size);
   ~Population();

   void rndInitialize(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
         vector<pair<double,double> > &range_realvar,
         vector<pair<double,double> > &range_bivar, Randomizer *r);

   void decode_pop();
   void evaluate_pop(MOP *mop);
   int getSize();
//   void copyPopulation(individual *source, int size);
};

#endif /* POBLACION_H_ */
