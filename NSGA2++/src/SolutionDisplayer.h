/*
 * SolutionDisplayer.h
 *
 *  Created on: 11/11/2009
 *      Author: antonio
 */

#ifndef SOLUTIONDISPLAYER_H_
#define SOLUTIONDISPLAYER_H_

#include <cstdio>
#include <fstream>
#include "Population.h"
#include "Individual.h"

class SolutionDisplayer {
public:
   void report_pop(Population *pop, ostream &output);
   void report_feasible (Population *pop, ostream &output);
   void report_ind(Individual *ind, ostream &output);
};

#endif /* SOLUTIONDISPLAYER_H_ */
