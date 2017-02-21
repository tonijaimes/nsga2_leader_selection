#ifndef SENSOR_H_
#define SENSOR_H_

#include "mop.h"
#include <vector>

class SensorNet : public MOP
{
public:
   SensorNet(int nNodos = 20, int nRows = 200, int nCols = 500);
	virtual ~SensorNet();
	
   void evaluate(vector<double> const &x, vector<int> const &gene,
                 vector<double> &fx, vector<double> &gcons) const;

//   void evaluate(vector<double> const &x, vector<double> &fx, vector<double> &gcons) const;
//   void evaluate(double const *x, double *fx, double *gcons) const;
//
//   void evaluate(vector<int> const &x, vector<double> &fx, vector<double> &gcons) const;

private:
   int const BASE_ID = 0;
   vector<pair<int,int> > sensorsPos;
   double const initEnergy = 9.0e6; // Energy given in micro Joules, 9x10^6 micro Joules = 9 Joules
   int numPackets = 10;


   double dist(pair<int,int> const &p1, pair<int,int> const &p2) const;
   double distPow(pair<int, int> const &x, pair<int, int> const &y) const;
   int findClosestLeader(vector<vector<int> > clusters, int i) const;
};

#endif /*SENSOR_H_*/
