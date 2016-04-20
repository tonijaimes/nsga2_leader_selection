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

using namespace std;


class NSGA {

public:
   //FILE *fpt1, *fpt2, *fpt3, *fpt4, *fpt5;
   ofstream file_firstPop;    //fpt1
   ofstream file_everyPop;    //fpt4
   ofstream file_finalPop;    //fpt2
   ofstream file_feasiblePop; //fpt3
   ofstream file_params;      //ftp5
   string name_firstPop;
   string name_everyPop;
   string name_finalPop;
   string name_feasiblePop;
   string name_params;

   MOP *mop; // Reference to MOP to be solve.
   Randomizer *r; // Reference to the random engine to be used.

   /* Attributes of NSGA */
   Population *parent_pop;
   Population *child_pop;
   Population *mixed_pop;

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


   NSGA(MOP *mop, Randomizer *r,
        int popsize, int ngen,
        int nreal, double pcross_real, double pmut_real, double eta_c, double eta_m,
        int nbin, double pcross_bin, double pmut_bin, vector<int> &nbits);
   ~NSGA();

   void optimize();

   void nGenerations(int n);

   void allocate_memory();
   void allocate_memory_pop (Population *pop, int size);
   void allocate_memory_ind (Individual *ind);
   void deallocate_memory_pop (Population *pop, int size);
   void deallocate_memory_ind (Individual *ind);

   double maximum (double a, double b);
   double minimum (double a, double b);

   void crossover (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
   void realcross (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);
   void bincross (Individual *parent1, Individual *parent2, Individual *child1, Individual *child2);

   void assign_crowding_distance_list (Population *pop, list *lst, int front_size);
   void assign_crowding_distance_indices (Population *pop, int c1, int c2);
   void assign_crowding_distance (Population *pop, int *dist, int **obj_array, int front_size);

   void decode_pop (Population *pop);
   void decode_ind (Individual *ind);

   int check_dominance (Individual *a, Individual *b);

   void evaluate_pop (Population *pop);
   void evaluate_ind (Individual *ind);

   void fill_nondominated_sort (Population *mixed_pop, Population *new_pop);
   void crowding_fill (Population *mixed_pop, Population *new_pop, int count, int front_size, list *cur);

   void initialize_pop (Population *pop);
   void initialize_ind (Individual *ind);

   void insert (list *node, int x);
   list* del (list *node);

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
   void quicksort_dist(Population *pop, int *dist, int front_size);
   void q_sort_dist(Population *pop, int *dist, int left, int right);

   void selection (Population *old_pop, Population *new_pop);
   Individual* tournament (Individual *ind1, Individual *ind2);

   void write_feasible(ostream &output) const;
   void write_whole_pop(ostream &output) const;

   void report_feasible(Population *pop, ostream &output) const;
   void report_pop(Population *pop, ostream &output) const;
   void report_params(ostream &output) const;
};

#endif /* NSGA_H_ */
