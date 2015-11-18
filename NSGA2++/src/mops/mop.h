#ifndef MOP_H
#define MOP_H

#include <string>
#include <vector>
#include <utility>

using namespace std;

/**
 * Abstract base class MOP defines a general Multiobjective Optimization Problem (MOP).
 * The abstract class MOP should be implemented (by inheritance) by any class
 * intended to define a concrete MOP. The child class must define the
 * evaluate() method, and, if needed, the evalConstraints() method .
 */
//template <class T>
class MOP {
public:
   MOP(const char *name, int numVariables, int numObjectives, int numConstraints = 0);
   virtual ~MOP();
   
   /** The next 2 pure virtual methods must be implemented in the actual MOP class */
   virtual void evaluate(vector<double> const &x, vector<double> &fx, vector<double> &gcons) const = 0;
   virtual void evaluate(double const *x, double *eval, double *gcons) const = 0;
   /********************************************************************************/

   virtual void setNumObjectives(int numObjectives);
   virtual void setNumVariables(int numVariables);
   virtual void setVarsRange(vector<pair<double,double> > &xRanges);
   //virtual void setEqualVarsRange(double lower, double upper);

   string getName() const;
   const int getNumVariables() const;
   const int getNumObjectives() const;
   const int getNumConstraints() const;
   const bool isScalableObjSpace() const;
   const bool isScalableVarSpace() const;
   vector<pair<double,double> > getVarsRange() const;


protected:
   string name;           /**< The name given to the MOP for information purposes. */
   int nVariables;        /**< Number of decision variables of the MOP. */
   int nObjectives;       /**< Number of objectives of the MOP. */
   int nConstraints;      /**< Number of constraints of the MOP. */
   bool isScalableObjSpc; /**< Indicates if the MOP is scalable in the objective space. */
   bool isScalableVarSpc; /**< Indicates if the MOP is scalable in the decision space. */
   vector<pair<double,double> > xRanges; /**< Range of the variables of the MOP. */

//   vector<double> lowerB; /**< Lower bounds for each variable of the MOP. */
//   vector<double> upperB; /**< Upper bounds for each variable of the MOP. */
}; // end abstract class MOP

#endif
