/*
 * VectorTools.h
 *
 *  Created on: 06-ago-2009
 *      Author: antonio
 */

#ifndef VECTORTOOLS_H_
#define VECTORTOOLS_H_


#include <fstream>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <utility>
#include <string>
#include <typeinfo>
#include <valarray>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

//namespace flab {

class VectorTools {
public:
	template <typename T>
   static string vectorToStr(vector<T> v,
                             int digitsPrec = 5,
                             const char *delimiter = ", ",
                             const pair<string, string> &enclosers = make_pair("[", "]")
                            );

	template <typename T>
   static void printVector(vector<T> v);

	template <typename T>
   static string matrixToStr(vector<vector<T> > m, int digitsPrec = 5, const char *delimiter = ", ");

	static vector<valarray<double> > loadDelimArray(string fileName, string delims ="\\s\\t");
	static vector<vector<double> >  loadDelimVector(string fileName, string delims ="\\s\\t");
};

/////////////////////////////////////////////////
// Implementation of the methods.
template <typename T>
void VectorTools::printVector(vector<T> v)
{
   cout << vectorToStr(v, 5, ", ");
}

template <typename T>
string VectorTools::vectorToStr(vector<T> v,
		              int digitsPrec,
		              const char *delimiter,
		              pair<string, string> const  &enclosers)
{
   if (v.size() == 0)
	   return "";

   if (digitsPrec < 0)
   	digitsPrec = 0;

	int gap;
   if (typeid(T) == typeid(int)) {
   	gap = 1;
   	digitsPrec = 0;
   }
   else
   	gap = 5;


   ostringstream strStream;
   strStream << fixed;
   strStream << enclosers.first;
   for (unsigned i = 0; i < v.size()-1; ++i) {
      strStream << setw(gap + digitsPrec) << setprecision(digitsPrec) << v[i] << delimiter;
   }

   strStream << fixed << setw(gap + digitsPrec) << setprecision(digitsPrec) << v.back();
   strStream << enclosers.second;

   return strStream.str();
}

template <typename T>
string VectorTools::matrixToStr(vector<vector<T> > m, int digitsPrec, const char *delimiter)
{
   ostringstream strStream;
   for (unsigned i = 0; i < m.size(); ++i)
      strStream << VectorTools::vectorToStr(m.at(i), digitsPrec, delimiter) << "\n";

   return strStream.str();
}

//}  // namespace flab

#endif /* VECTORTOOLS_H_ */
