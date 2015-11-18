/*
 * RandUtils.h
 *
 *  Created on: 30/10/2009
 *      Author: antonio
 */

#ifndef RANDUTILS_H_
#define RANDUTILS_H_


class Randomizer {
public:
   Randomizer();
   Randomizer(double seed);
   ~Randomizer();

   /* Function declarations for the random number generator */
    void setSeed(double seed);
    double getSeed();
    void randomize();
    void warmup_random(double seed);
    void advance_random(void);
    double randomperc(void);
    int rnd(int low, int high);
    double rndreal(double low, double high);

private:
   /* Variable declarations for the random number generator */
   double seed;
   double oldrand[55];
   int jrand;
};

#endif /* RANDUTILS_H_ */
