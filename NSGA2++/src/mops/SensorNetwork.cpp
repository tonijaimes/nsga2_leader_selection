#include "SensorNetwork.h"
#include <cmath>
#include <iostream>

using namespace std;


SensorNet::SensorNet(int nNodes, int nRows, int nCols) : MOP("SensorNet", nNodes, 2, 0)
{
   // The base is located in the corner (0,0) of the grid.
   sensorsPos.push_back(make_pair(0, 0));

   srand(100);

   // Generate random positions for nNodes sensors in a nRows x nCols grid
   int x, y, i;
   for (i = 0; i < nNodes; ++i) {
      x = 1 + rand()%nCols ;// Random integer in 1,...,nCols
      y = 1 + rand()%nRows; // Random integer in 1,...,nRows
      sensorsPos.push_back( make_pair(x, y) );
   }
}

SensorNet::~SensorNet() {}

/*
 * leaderOrNot: binary vector indicating whether i-th node is leader (value 1) or not (value 0)
 * eval: a vector in which each computed objective value will be stored.
 * gcons: to store constraint values, but it is not used here.
 */
void SensorNet::evaluate(vector<int> const &leaderOrNot, vector<double> &eval,  vector<double> &gcons) const
{

   // Gather information of the clusters formed.
   // Each cluster vector will contain the indexes of the nodes conforming the cluster.
   vector<vector<int> > clusters;
   vector<double> energy(leaderOrNot.size(), initEnergy); // Residual energy of each sensor.

   //First only the leaders of each cluster are inserted.
   int pos = 1; // positions start at 1, because the base is located at position 0.
   int numLeaders = 0;
   for (auto isLeader : leaderOrNot) {
      if (isLeader) {
         clusters.emplace_back( pos );
         numLeaders++;
      }
      pos++;
   }

   // Obtain the members of each cluster using k-medoids method.
   // I.e., which is the closest leader for each member.
   pos = 1;
   int myCluster;
   for (auto isLeader : leaderOrNot) {
      if (!isLeader) {
    	  myCluster = findClosestLeader(clusters, pos);
        clusters[myCluster].push_back(pos);
      }
      pos++;
   }


   // Objective 1: Distance of leaders to the base
   double totalDist = 0;
   pos = 1;
   for (auto c : clusters) {
      totalDist += dist(sensorsPos[ c[0] ], sensorsPos[0]);
   }
   eval[0] = totalDist / numLeaders;


   // Objective 2: Distance of cluster members to their leader.
   double sumDistCluster;
   double sumDist;
   for (auto c : clusters) {

      sumDistCluster = 0;
      for (unsigned int i = 1; i < c.size(); ++i) {
         // Sumar distancia del leader al sensor i
         sumDistCluster += dist(sensorsPos[ c[0] ], sensorsPos[ c[i] ]);
      }
      sumDist += sumDistCluster / c.size();
   }
   eval[1] = sumDist / numLeaders;



   // Objective 3: Residual energy.
   // Compute energy used for the current round of messaging
   // For leaders
   double energyUsed;
   for (auto c : clusters) {
      energyUsed = 100 +
                   0.2*distPow(sensorsPos[ c[0] ], sensorsPos[BASE_ID]) +
                   10 *(c.size()-1);
      energy[ c[0] ] -= energyUsed;

      //if (energy[ c[0] ] <= 0) {
      //   energy[ c[0] ] = 0.0;
      //}
   }


   // For cluster members
   for (auto c : clusters) {
      for (unsigned int i = 1; i < c.size(); ++i) {
         energyUsed = 100 + 0.2* distPow(sensorsPos[ c[i] ], sensorsPos[ c[0] ]);
         energy[ c[0] ] -= energyUsed;
      //   if (energy[i] <= 0) {
         //}
      }
   }

   double energyLeaders=0;
   for (auto c : clusters) {
      energyLeaders += energy[ c[0] ];
   }

   eval[2] = energyLeaders / numLeaders;

}

double SensorNet::dist(pair<int,int> const &x, pair<int,int> const &y) const {
   double d = sqrt( pow(x.first - y.first, 2) + pow(x.second - y.second, 2) );

   return d;
}

double SensorNet::distPow(pair<int, int> const &x, pair<int, int> const &y) const
{
   double d = pow(x.first - y.first, 2) + pow(x.second - y.second, 2);

   return d;
}

int SensorNet::findClosestLeader(vector<vector<int> > clusters, int i) const {
	int myCluster;

	int currentCluster = 0;
	double minDist = 1e7;
	double d;
	for (auto c : clusters) {
	   // Compute distance from leader of cluster c to i-th node.
	   d = dist(sensorsPos[ c[0] ], sensorsPos[i]);
	   if (d < minDist) {
	      minDist = d;
	      myCluster = currentCluster;
	   }
	   currentCluster++;
	}

	return myCluster;
}


void SensorNet::evaluate(double const *x, double *eval, double *gcons) const {}

void SensorNet::evaluate(vector<double> const &x, vector<double> &eval,  vector<double> &gcons) const {}

