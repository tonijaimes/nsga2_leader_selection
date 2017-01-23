/*
 * NSGA.h
 *
 *  Created on: 28/10/2009
 *      Author: antonio
 */

#ifndef NSGA_H_
#define NSGA_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "RandUtils.h"
#include "List.h"
#include "Population.h"
#include "Individual.h"
#include "theMOPS.h"
//#include "KOSSA.h"
//#include "Partition.h"
#include "myUtils.h"
///#include "Archive.h"

using namespace std;

/*
 * Identifies the type of search phase performed:
 * integrate all the objectives, or search on
 * separate objective subspaces.
 */
enum PhasesType
{
   INTEGRATION,
   PARTITIONING
};

/* Contains all the possible parameters related
 * to the partitioning method.
 */
struct PartitionParams {
	int numSubspaces;
	int delta0;
	int numPhases;
	double incrFactor;
	double pPhase;
	double pPartition;
	string partitionType;
};

/* Contains all the possible parameters of
 * NSGA.
 */

struct NSGAParams {
   int popsize;        // Population size.
   int ngen;           // Number of generations.
   int nreal;          // Number of variables using a real encoding.
   double pcross_real; // Crossover probability for real variables.
   double pmut_real;   // Mutation probability for real variables.
   double eta_c;
   double eta_m;
   int nbin;           // Number of variables using binary encoding.
   vector<int> nbits;  // Bits for each binary encoded variable.
   double pcross_bin;  // Crossover probability for binary variables.
   double pmut_bin;    // Mutation probability for binary variables.
};

class NSGA {

public:
	bool verbose;
   //FILE *fpt1, *fpt2, *fpt3, *fpt4, *fpt5;
   ofstream file_firstPop;    //fpt1
   ofstream file_everyPop;    //fpt4
   ofstream file_finalPop;    //fpt2
   ofstream file_feasiblePop; //fpt3
   ofstream file_params;      //ftp5

   ofstream file_varSpace;
   ofstream file_objSpace;
   ofstream file_everySubParentPop;
   ofstream file_everySubChildPop;
   ofstream file_inParentPop;
   ofstream file_subspaces;
   ofstream file_subspacesInfo;
   ofstream file_onlineInd;
   ofstream file_archive;

   string name_firstPop;
   string name_everyPop;
   string name_finalPop;
   string name_feasiblePop;
   string name_params;

   string name_everySubParentPop;
   string name_everySubChildPop;
   string name_varSpace;
   string name_objSpace;
   string name_inParentPop;
   string name_subspaces;
   string name_subspacesInfo;
   string name_onlineInd;

   //Archivos para el esquema adaptativo.
   ofstream file_partPhaseLength;
   ofstream file_intePhaseLength;
   ofstream file_convSubspaces;
   ofstream file_convWholeSpace;
   ofstream file_parentShare;

   MOP *mop;
   Randomizer *r;

   /* Atributos del NSGA */
   Population *parent_pop;
   Population *child_pop;
   Population *mixed_pop;

//   Archive *archive;

   int ngen;
   int popsize;
   int nbinmut;
   int nrealmut;
   int nbincross;
   int nrealcross;
   vector<int> nbits;
   vector<pair<double,double> > range_realvar;
   vector<pair<double,double> > range_binvar;
   int bitlength;


   /* Atributos del MOP */
   int nreal;
   int nbin;
   int nobj;
   int ncon;

   /* Parámetros operadores */
   double pcross_real;
   double pcross_bin;
   double pmut_real;
   double pmut_bin;
   double eta_c;
   double eta_m;

   PartitionParams partition;
   vector<int> entireObjSet;
//   vector<doubleRange >  range_ndset;
   bool rndPopulation;
   vector<vector<double> > initialPop;


   NSGA(MOP *mop, Randomizer *r, PartitionParams &p, NSGAParams &np);

   ~NSGA();

   void startWithThisPopulation(vector<vector<double> > const &initPop);
   void startWithRandomPopulation();

   void archivePopulation(Population *pop);
//   void initPFRange(Population *pop);
//   void updatePFRange(Population *pop);
   string rangesToStr(int digitsPrec = 5, const char *delimiter = ", ");
//   void normalizeNDSet(doubleMatrix &ndsets, vector<int> &objSet, vector<doubleRange > &range_ndset);
//   void normalizeMultiNDSet(vector<doubleMatrix > &ndsets, intMatrix &objSpaces, vector<doubleRange > &range_ndset);
//   vector<doubleMatrix > recordNonDominatedSets(Population *pop, intMatrix partition);
//	void updatePartition(Partition &partition, vector<bool> convergenceResults);
//	bool allSubspacesConverged(vector<bool> hasSubspaceConverged);
//	bool someSubspaceConverged(vector<bool> hasConverged);

   void createInitialPopulation();
   void optimize();

   void openInfoFiles();
   void closeInfoFiles();

   void nGenerations(int n);
   /*
   vector<vector<int> > createPartition(Population *parent_pop, ClusterInfo &info);
   Partition createPartition(Population *parent_pop);
   void subspaces_fill_nondominated_sort(Population *mixed_pop, Population *new_pop, vector<vector<int> > &oSpaces);
   void proportional_subspaces_fill_nondominated_sort(
   		Population *mixed_pop, Population *new_pop, Partition &partition);
   void proportional_subspaces_fill_nondominated_sort(Population *mixed_pop, Population *new_pop,
   		                                             vector<vector<int> > &oSpaces, ClusterInfo &info);
   void subspaces_SortAndTruncationAndSelection(Population *mixed_pop, Population *new_pop, Population *child_pop, vector<vector<int> > &oSpaces);
   vector<vector<int> > setOSpaces_Fixed(int numObjs, int numSubSpaces);
   vector<vector<int> > setOSpaces_Randomly(int numObjs, int numSubSpaces);
   void displaySubSpaces(vector<vector<int> > &oSpaces, ostream &output,
                         char const *left="[", char const *right="]", char const *separator=" ");
   void displaySubSpacesInfo(vector<vector<int> > &oSpaces, ClusterInfo &info);
   void reportSubSpacesInfo(vector<vector<int> > &oSpaces, ClusterInfo &info, ostream &output);
*/
   void allocate_memory();
   void allocate_memory_pop (Population *pop, int size);
   void allocate_memory_ind (Individual *ind);
   void deallocate_memory_pop (Population *pop, int size);
   void deallocate_memory_ind (Individual *ind);

   void crossover (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
   void realcross (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
   void bincross (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);

   void assign_crowding_distance_list(Population *pop, myList *lst, int front_size);
   void assign_crowding_distance_list(Population *pop, myList *lst, int front_size, vector<int> &objSet);
   void assign_crowding_distance_indices (Population *pop, int c1, int c2);
   void assign_crowding_distance_indices (Population *pop, int c1, int c2, vector<int> &objSet);
   void assign_crowding_distance (Population *pop, vector<int> &dist, vector<vector<int> >  &obj_array,
                                  vector<int> &objSet, int front_size);
   void assign_crowding_distance (Population *pop, int *dist, int **obj_array, int front_size);

   void decode_pop (Population *pop);
   void decode_ind (Individual *ind);

   int check_dominance (Individual *a, Individual *b);
   int check_dominance (Individual *a, Individual *b, vector<int> &objSet);

   void evaluate_pop (Population *pop);
   void evaluate_ind (Individual *ind);

   void fill_nondominated_sort (Population *mixed_pop, Population *new_pop, vector<int> &objSet);
   void fill_nondominated_sort (Population *mixed_pop, Population *new_pop);
   void crowding_fill(Population *mixed_pop, Population *new_pop, int count, int front_size,
                      myList *cur, vector<int> &objSet);
   void crowding_fill(Population *mixed_pop, Population *new_pop, int count, int front_size,
   		             myList *elite);

   void initialize_pop (Population *pop);
   void initialize_pop_fix(Population *pop);
   void initialize_ind (Individual *ind);

   void insert (myList *node, int x);
   myList* del (myList *node);

   void merge(Population *pop1, Population *pop2, Population *pop3);
   void copy_ind (Individual *ind1, Individual *ind2);

   void mutation_pop (Population *pop);
   void mutation_ind (Individual *ind);
   void bin_mutate_ind (Individual *ind);
   void real_mutate_ind (Individual *ind);

   void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr);
   void assign_rank_and_crowding_distance (Population *new_pop);

   void quicksort_front_obj(Population *pop, int objcount, int obj_array[], int obj_array_size);
   void q_sort_front_obj(Population *pop, int objcount, int obj_array[], int left, int right);
   void quicksort_front_obj(Population *pop, int objcount, vector<int> &obj_array, int obj_array_size);
   void q_sort_front_obj(Population *pop, int objcount, vector<int> &obj_array, int left, int right);

   void quicksort_dist(Population *pop, vector<int> &dist, int front_size);
   void q_sort_dist(Population *pop, vector<int> &dist, int left, int right);
   void quicksort_dist(Population *pop, int *dist, int front_size);
   void q_sort_dist(Population *pop, int *dist, int left, int right);

   void selection (Population *old_pop, Population *new_pop);
   Individual* tournament (Individual *ind1, Individual *ind2);

   void write_feasible(ostream &output) const;
   void write_whole_pop(ostream &output) const;

   void report_feasible(Population *pop, ostream &output) const;
   void report_objSpace(Population *pop, ostream &output) const;
   void report_varSpace(Population *pop, ostream &output) const;
   void report_pop(Population *pop, ostream &output) const;
   void report_params(ostream &output) const;
   void report_pop_objs(Population *pop, ostream &output) const;
};

#endif /* NSGA_H_ */
