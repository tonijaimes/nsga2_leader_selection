/*
 * RandUtils.cpp
 *
 *  Created on: 30/10/2009
 *      Author: antonio
 */

#include "RandUtils.h"
#include <ctime>
#include <cstdlib>
#include <cassert>

using namespace std;

Randomizer::Randomizer() {
   srand( time(0) );
   setSeed( ((double) random()) / RAND_MAX );
}

Randomizer::Randomizer(double seed) {
   setSeed(seed);
}

Randomizer::~Randomizer() {}

void Randomizer::setSeed(double seed) {
   assert(seed <= 0 || seed >= 1.0);

   this->seed = seed;
}

double Randomizer::getSeed() {
   return seed;
}

/* Get seed number for random and start it up */
void Randomizer::randomize()
{
    int j1;
    for(j1=0; j1<=54; j1++)
    {
        oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random(seed);
    return;
}

/* Get randomize off and running */
void Randomizer::warmup_random (double seed)
{
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = seed;
    new_random = 0.000000001;
    prev_random = seed;
    for(j1=1; j1<=54; j1++)
    {
        ii = (21*j1)%54;
        oldrand[ii] = new_random;
        new_random = prev_random-new_random;
        if(new_random<0.0)
        {
            new_random += 1.0;
        }
        prev_random = oldrand[ii];
    }
    advance_random ();
    advance_random ();
    advance_random ();
    jrand = 0;
    return;
}

/* Create next batch of 55 random numbers */
void Randomizer::advance_random ()
{
    int j1;
    double new_random;
    for(j1=0; j1<24; j1++)
    {
        new_random = oldrand[j1]-oldrand[j1+31];
        if(new_random<0.0)
        {
            new_random = new_random+1.0;
        }
        oldrand[j1] = new_random;
    }
    for(j1=24; j1<55; j1++)
    {
        new_random = oldrand[j1]-oldrand[j1-24];
        if(new_random<0.0)
        {
            new_random = new_random+1.0;
        }
        oldrand[j1] = new_random;
    }
}

/* Fetch a single random number between 0.0 and 1.0 */
double Randomizer::randomperc()
{
    jrand++;
    if(jrand>=55)
    {
        jrand = 1;
        advance_random();
    }
    return((double)oldrand[jrand]);
}

/* Fetch a single random integer between low and high including the bounds */
int Randomizer::rnd(int low, int high)
{
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        res = low + (randomperc()*(high-low+1));
        if (res > high)
        {
            res = high;
        }
    }
    return (res);
}

/* Fetch a single random real number between low and high including the bounds */
double Randomizer::rndreal(double low, double high)
{
    return (low + (high-low)*randomperc());
}
