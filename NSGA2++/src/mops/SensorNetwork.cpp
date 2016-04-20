#include "SensorNetwork.h"
#include <cmath>
#include <iostream>

using namespace std;


SensorNet::SensorNet(int nNodes, int nRows, int nCols) : MOP("SensorNet", nNodes, 2, 0)
{
   // The base is located in the corner (0,0) of the grid.
   sensorsPos.push_back(make_pair(0, 0))

   srand();

   // Generate random positions for nNodes sensors in a nRowsxnCols grid
   int x, y;
   for (i = 0; i < nNodes; ++i) {
      x = 1 + rand()%nCols ;// Random integer in 1,...,nCols
      y = 1 + rand()%nRows; // Random integer in 1,...,nRows
      sensorsPos.push_back( make_pair(x, y) );
   }
}

SensorNet::~SensorNet() {}

void SensorNet::evaluate(vector<int> const &leaderOrNot, vector<double> &eval,  vector<double> &gcons) const
{
   // Gather informatio of the clusters formed.
   // Obtain the leaders encoded in the binaryString
   int pos = 1;
   int numLeaders = 0;
   for (isLeader : leaderOrNot) {
      if (isLeader) {
         clusters.push_back( vector<int>(pos) );
         numLeaders++;
      }
      pos++;
   }

   // Obtain the members of each cluster using k-medoids
   // I.e., which is the closest leader for each member.
   pos = 1;
   for (isLeader : leaderOrNot) {
      if (!isLeader) {
         myLeader = findClosestLeader(clusters);
         clusters[myLeader].push_back(pos);
      }
      pos++;
   }


   // Objective 1: Distance of leaders to the base
   double totalDist = 0;
   int pos = 1;
   for (c : clusters) {
      totalDist += dist(sensorsPos[ c[0] ], sensorsPos[0]);
   }
   eval[0] = totalDist / numLeaders;


   // Objective 2: Distance of cluster members to their leader.
   double sumDistCluster;
   double sumDist;
   for (c : clusters) {

      sumDistCluster = 0;
      for (int i = 1; i < c.size(); ++i) {
         // Sumar distancia del leader al sensor i
         sumDistCluster += dist(sensorsPos[ c[0] ], sensorsPos[ [i] ]);
      }
      sumDist += sumDistCluster;
   }
   eval[1] = sumDist / numLeaders;



   // Objective 3: Residual energy.
   eval[2] = 0;

}

void SensorNet::evaluate(double const *x, double *eval, double *gcons) const {}

void SensorNet::evaluate(vector<double> const &x, vector<double> &eval,  vector<double> &gcons) const {}

