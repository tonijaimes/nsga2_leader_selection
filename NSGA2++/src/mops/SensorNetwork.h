#ifndef SENSOR_H_
#define SENSOR_H_

#include "mop.h"
#include <vector>

class SensorNet : public MOP
{
public:
   SensorNet(int nNodos = 20, int nRows = 20, int nCols = 50);
	virtual ~SensorNet();
	
   void evaluate(vector<double> const &x, vector<double> &fx, vector<double> &gcons) const;
   void evaluate(double const *x, double *fx, double *gcons) const;

   void evaluate(vector<int> const &x, vector<double> &fx, vector<double> &gcons) const;

private:
   vector<pair<int,int> > sensorsPos;
   double dist(pair<int,int> const &p1, pair<int,int> const &p2) const;
   int findClosestLeader(vector<vector<int> > clusters, int i) const;
};

#endif /*SENSOR_H_*/