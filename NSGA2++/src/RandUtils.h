/*
 * RandUtils.h
 *
 *  Created on: 30/10/2009
 *      Author: antonio
 */

#ifndef RANDUTILS_H_
#define RANDUTILS_H_

#include <vector>

using namespace std;

class Randomizer {
public:
   Randomizer();
   Randomizer(double seed);
   ~Randomizer();

   /* Function declarations for the random number generator */
    void setSeed(double seed);
    double getSeed();
    void randomize();
    double randomperc(void);
    int rnd(int low, int high);
    double rndreal(double low, double high);
    void rndVector(vector<int> &rndVector);

private:
   static const int oldRandSize; /* One for all the instances. */
   /* Variable declarations for the random number generator */
   double seed;
   vector<double> oldrand;
   int jrand;

   //double oldrand[55];

   void warmup_random(double seed);
   void advance_random(void);
};

#endif /* RANDUTILS_H_ */
