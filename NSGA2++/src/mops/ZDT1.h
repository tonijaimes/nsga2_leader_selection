#ifndef ZDT1_H_
#define ZDT1_H_

#include "mop.h"
#include <vector>

class ZDT1 : public MOP
{
public:
	ZDT1(int nVariables = 30);
	virtual ~ZDT1();
	
   void evaluate(vector<double> const &x, vector<double> &fx, vector<double> &gcons) const;
   void evaluate(double const *x, double *fx, double *gcons) const;
};

#endif /*ZDT1_H_*/
