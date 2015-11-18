#include <iostream>
//#include <mop.h>
#include "mop.h"

using std::string;

MOP::MOP(const char *name, int numVariables, int numObjectives, int numConstraints)
{
   this->name.assign(name);
   this->nObjectives = numObjectives;
   this->nVariables = numVariables;
   this->nConstraints = numConstraints;
   isScalableObjSpc = false;
   isScalableVarSpc = false;

   // These are some default values, each concrete MOP should set appropriate values.
   xRanges.assign(numVariables, make_pair(0.0, 1.0));
}

MOP::~MOP() {}

void MOP::setNumObjectives(int numObjectives) {
   if (isScalableObjSpc)
      this->nObjectives = numObjectives;
}

void MOP::setNumVariables(int numVariables) {
   if (isScalableVarSpc)
      this->nVariables = numVariables;
}

void MOP::setVarsRange(vector<pair<double,double> > &xRanges) {
   this->xRanges = xRanges;
}

const bool MOP::isScalableObjSpace() const{
   return isScalableObjSpc;
}

const bool MOP::isScalableVarSpace() const{
   return isScalableVarSpc;
}

const int MOP::getNumVariables() const {
   return nVariables;
}

vector<pair<double,double> > MOP::getVarsRange() const {
   return xRanges;
}

const int MOP::getNumObjectives() const {
   return nObjectives;
}

const int MOP::getNumConstraints() const {
   return nConstraints;
}

string MOP::getName() const {
   return name;
}
