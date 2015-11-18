/*
 * Individual.h
 *
 *  Created on: 29/10/2009
 *      Author: antonio
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

/*
typedef struct
{
    int rank;
    double constr_violation;
    double *xreal;
    int **gene;
    double *xbin;
    double *obj;
    double *constr;
    double crowd_dist;

    int numObjs;
    int numVars;
    int numCons;
}
individual;
*/

#include <vector>
#include "RandUtils.h"


using namespace std;

class Individual {
public:
   Individual();
   Individual(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
         vector<pair<double,double> > &range_realvar,
         vector<pair<double,double> > &range_bivar, Randomizer *r);
   ~Individual();

   void rndInitialize(int nreal, int nobj, int ncon, int nbin, vector<int> &nbits,
         vector<pair<double,double> > &range_realvar,
         vector<pair<double,double> > &range_bivar, Randomizer *r);

   int rank;
   double constr_violation;
   double crowd_dist;
   vector<double> obj;
   vector<double> constr;
   vector<double> xreal;
   vector<double> xbin;
   int nbin;
   vector<vector<int> > gene;
   vector<int> nbits;
   vector<pair<double,double> >  range_realvar;
   vector<pair<double,double> >  range_binvar;

   void decode();

   void setRank(double rank);
   int getRank();

   void setCrowd_dist(double crowd_dist);
   double getCrowd_dist();

   void setConstr_violation(double constr_viol);
   double getConstr_violation();
};

#endif /* INDIVIDUAL_H_ */
