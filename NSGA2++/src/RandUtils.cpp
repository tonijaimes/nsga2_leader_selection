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
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>

using namespace std;

const int Randomizer::oldRandSize = 55;

Randomizer::Randomizer() {
   oldrand.resize(oldRandSize, 0.0);
   srand( time(0) );
   setSeed(((double) random()) / RAND_MAX);
}

Randomizer::Randomizer(double seed) {
   oldrand.resize(oldRandSize, 0.0);
   setSeed(seed);
}

Randomizer::~Randomizer() {}

void Randomizer::setSeed(double seed) {
   assert(seed > 0.0 && seed < 1.0);
   this->seed = seed;
}

double Randomizer::getSeed() {
   return seed;
}

/* Get seed number for random and start it up */
void Randomizer::randomize()
{
   oldrand.assign(oldRandSize, 0.0);
   jrand = 0;
   warmup_random(seed);
}

/* Get randomize off and running */
void Randomizer::warmup_random(double seed)
{
   unsigned j1, ii;
   double new_random, prev_random;

   oldrand.back() = seed;
   //    oldrand[54] = seed;
   new_random = 0.000000001;
   prev_random = seed;
   for(j1=1; j1 < oldrand.size(); ++j1)
   {
      ii = (21*j1) % 54;
      //cout << ii << " ";
      try {
         oldrand.at(ii) = new_random;
         //oldrand[ii] = new_random;
         new_random = prev_random - new_random;
         if(new_random < 0.0)
            new_random += 1.0;

         prev_random = oldrand.at(ii);
      }
      catch ( std::out_of_range outOfRange ) {// out_of_range exception
         cout << "\n\nException: " << outOfRange.what();
      }
   }
   advance_random ();
   advance_random ();
   advance_random ();
   jrand = 0;

   /*
    cout << "\n" << endl;
    int jj;
    cout << fixed << setprecision(8);
    for(jj=0; jj <= 54; ++jj)
       cout << oldrand[jj] << endl;
    */
}

/* Create next batch of 55 random numbers */
void Randomizer::advance_random ()
{
    int j1;
    double new_random;

    try {
       for(j1=0; j1 < 24; ++j1)
       {
          new_random = oldrand.at(j1) - oldrand.at(j1+31);
          //new_random = oldrand[j1] - oldrand[j1+31];
          if(new_random < 0.0)
             new_random = new_random + 1.0;

          oldrand.at(j1) = new_random;
          //oldrand[j1] = new_random;
       }
       for(j1=24; j1 < 55; ++j1)
       {
          new_random = oldrand.at(j1) - oldrand.at(j1-24);
          //new_random = oldrand[j1] - oldrand[j1-24];
          if(new_random<0.0)
             new_random = new_random + 1.0;

          oldrand.at(j1) = new_random;
          //oldrand[j1] = new_random;
       }
    }
    catch ( std::out_of_range outOfRange ) {// out_of_range exception
       cout << "\n\nException: " << outOfRange.what();
    }
}

/* Fetch a single random number between 0.0 and 1.0 */
double Randomizer::randomperc()
{
    jrand++;
    if(jrand >= 55)
    {
        jrand = 1;
        advance_random();
    }

    return((double) oldrand.at(jrand));
}

/* Fetch a single random integer between low and high including the bounds */
int Randomizer::rnd(int low, int high)
{
    int res;
    if (low >= high)
        res = low;
    else
    {
        res = low + (randomperc() * (high-low+1));
        if (res > high)
            res = high;
    }

    return res;
}

/* Fetch a single random real number between low and high including the bounds */
double Randomizer::rndreal(double low, double high)
{
    double perc = randomperc();
//    cout << fixed << setprecision(8);
//    cout << perc << endl;
    return (low + (high-low)*perc);
//    return (low + (high-low)*randomperc());
}

void Randomizer::rndVector(vector<int> &rndVec) {
   int rndPick;

   for (unsigned i = 0; i < rndVec.size(); ++i)
      rndVec.at(i) = i;

   unsigned last = rndVec.size()-1; /* But last index value */
   for (unsigned first = 0; first < last; ++first) {
      rndPick = rnd(first, last);
      swap(rndVec[rndPick], rndVec[first]);
   }
}

