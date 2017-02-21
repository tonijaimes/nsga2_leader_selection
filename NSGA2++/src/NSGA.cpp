/*
 * NSGA.cpp
 *
 *  Created on: 28/10/2009
 *      Author: antonio
 */
#include "Konstants.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include "NSGA.h"
//#include "MetricsProxy.h"
//#include "Metrics.h"
//#include "OStasisDetector.h"
//#include "VectorTools.h"
//#include "Archive.h"

using namespace std;

//#define CONFLICT


NSGA::NSGA(MOP *mop, Randomizer *r, NSGAParams &gaParams)
{
   this->mop = mop;
   this->r = r;
   this->popsize = gaParams.popsize;
   this->ngen = gaParams.ngen;
   this->nreal = gaParams.nreal;
   this->pcross_real = gaParams.pcross_real;
   this->pmut_real = gaParams.pmut_real;
   this->eta_c = gaParams.eta_c;
   this->eta_m = gaParams.eta_m;
   this->nbin = gaParams.nbin;
   this->pcross_bin = gaParams.pcross_bin;
   this->pmut_bin = gaParams.pmut_bin;
   this->nbits = gaParams.nbits;

   //partition = partitionParams;
   //archive = new Archive(1*popsize);

   /* The smallest subspaces must have at least size 2. */
   //assert(mop->getNumObjectives() / partition.numSubspaces >= 2);

   //cerr << "\nSize VarsRange= " << mop->getVarsRange().size() << endl;

   /* Take the variable range information from the MOP. */
   cout << "nreal = " << nreal << "\n";


   vector<pair<double,double> >::iterator inicio_realvar = (mop->xRanges).begin();
   vector<pair<double,double> >::iterator fin_realvar = (mop->xRanges).begin() + nreal;
   vector<pair<double,double> >::iterator fin_binvar = (mop->xRanges).end();

   if (nreal != 0) {
   	range_realvar.assign(inicio_realvar, fin_realvar);
   	cerr << "\nSize range_realvar= " << range_realvar.size() << endl;
   	cerr << "first interval: [" << range_realvar[0].first << ", " << range_realvar[0].second << "].\n";
   }

   if (nbin != 0) {
   	range_binvar.assign(fin_realvar, fin_binvar);
   	cerr << "\nSize range_binvar= " << range_binvar.size() << endl;
   }

   nobj = mop->getNumObjectives();
   ncon = mop->getNumConstraints();

   /* The entire objective space, Phi = {1,...,nobj}. */
   for (int i = 0; i < nobj; ++i)
   	entireObjSet.push_back(i);

   bitlength = 0;
   if (nbin != 0) {
       for (int i=0; i < nbin; ++i)
          bitlength += nbits[i];
   }

   nbinmut = 0;
   nrealmut = 0;
   nbincross = 0;
   nrealcross = 0;

   parent_pop = new Population(popsize);
   child_pop  = new Population(popsize);
   mixed_pop  = new Population(2*popsize);

   name_firstPop = "first_pop.out";
   name_everyPop = "all_pop.out";
   name_finalPop = "final_pop.out";
   name_feasiblePop = "best_pop.out";
   name_params = "params.out";

   //name_everySubParentPop = "all_subparents.out";
   //name_everySubChildPop = "all_subchildren.out";
   name_varSpace = "var_space.out";
   name_objSpace = "obj_space.out";
   name_inParentPop = "all_inpop.out";
   //name_subspaces = "subspaces.out";
   //name_subspacesInfo = "subInfo.out";
   name_onlineInd = "onlineIndicator.out";
}

void NSGA::allocate_memory() {
   /*
   allocate_memory_pop (parent_pop, popsize);
   allocate_memory_pop (child_pop, popsize);
   allocate_memory_pop (mixed_pop, 2*popsize);
   */
}

NSGA::~NSGA() {
	delete parent_pop;
	delete child_pop;
	delete mixed_pop;
}

void NSGA::startWithRandomPopulation() {
   rndPopulation = true;
}

void NSGA::startWithThisPopulation(vector<vector<double> > const &initPop) {
   cerr << "\nNumber of given solutions: " << initPop.size() << endl;

   assert((int)initPop.size() >= popsize);
   cerr << "\nA given pop will be used.\n" << endl;
   rndPopulation = false;
   this->initialPop = initPop;
}

void NSGA::nGenerations(int nGens) {

   for (int i=1; i <= nGens; i++)
   {
       selection (parent_pop, child_pop);
       mutation_pop (child_pop);
       decode_pop(child_pop);
       evaluate_pop(child_pop);
       merge (parent_pop, child_pop, mixed_pop);
       fill_nondominated_sort(mixed_pop, parent_pop);
       /* Comment following three lines if information for all
       generations is not desired, it will speed up the execution */
//       fprintf(fpt4,"# gen = %d\n",i);
//       report_pop(parent_pop,fpt4);
//       fflush(fpt4);
       printf("\n gen = %d",i);
   }
}

//void NSGA::archivePopulation(Population *pop) {
//	for (int i = 0; i < pop->getSize(); ++i)
//		archive->filter(*(pop->ind[i]));
//}

void NSGA::createInitialPopulation() {
   if (rndPopulation)
      initialize_pop(parent_pop);
   else
      initialize_pop_fix(parent_pop);

   initialize_pop(child_pop);
   initialize_pop(mixed_pop);
}

void NSGA::initialize_pop(Population *pop) {
   for (int i=0; i < pop->size; ++i)
      pop->ind[i]->rndInitialize(nreal, nobj, ncon, nbin, nbits, range_realvar, range_binvar, r);
//       initialize_ind (&(pop->ind[i]));
}

void NSGA::initialize_pop_fix(Population *pop) {
   for (int i=0; i < pop->size; ++i)
      pop->ind[i]->initialize(initialPop.at(i), nreal, nobj, ncon, nbin, nbits, range_realvar, range_binvar);
}

void NSGA::optimize()
{
	openInfoFiles();

   r->randomize();
   createInitialPopulation();

//   initialize_pop(parent_pop);
//   initialize_pop(child_pop);
//   initialize_pop(mixed_pop);

   cout << "\n Initialization done, now performing first generation.";

   //TODO: Para sensores no se necesita decodificar
   //decode_pop(parent_pop);

   evaluate_pop (parent_pop);
   assign_rank_and_crowding_distance(parent_pop);
   report_pop(parent_pop, file_firstPop);
   //initPFRange(parent_pop);
   file_everyPop << "# gen = 1\n";
   report_pop_objs(parent_pop, file_everyPop);
   file_everyPop << "# end gen = 1\n";
   cout << "\n gen = 1";

   int printingRate = 1 + ngen/10;
   //int maxSampleSize = 10; //TODO: this should be an input parameter.


//   Partition basicPartition(intMatrix(1, entireObjSet));
//   Partition currentPartition = basicPartition;
//   vector<doubleMatrix > oldNDSs, newNDSs;
//   vector<double> metricValues;
//   vector<bool> convergenceFlags;
//   OStasisDetector stasisDetector(currentPartition.oSpaces.size(), maxSampleSize, alphaPhi, sigLevel);
////   OStasisDetector stoppingDetector(basicPartition.oSpaces.size(), maxSampleSize, stopping, sigLevel);
//   OStasisDetector stoppingDetector(partition.numSubspaces, maxSampleSize, stopping, sigLevel);
//   //OStasisDetector stasisDector(sampleSize, 0.08, 0.05);
//   int currentNumSamples = 0; //Number of samples recorded.
//   PhasesType currentPhase = INTEGRATION;
//   int phaseGen = 1; //Generations used in the current phase (Integration or Partitioning).
   int i;
   verbose = false;

//   archivePopulation(parent_pop);

   for (i=2; i <= ngen; ++i)
   {
      // Dominance-based selection and children generation by crossover.
      selection(parent_pop, child_pop);

      mutation_pop(child_pop);
      decode_pop(child_pop);
      evaluate_pop(child_pop);
      merge(parent_pop, child_pop, mixed_pop); // Copy child an parent pop into mixed pop.

      fill_nondominated_sort(mixed_pop, parent_pop);

      file_everyPop << "# gen = " << i << "\n";
      report_pop_objs(parent_pop, file_everyPop);
      file_everyPop << "# end gen = " << i << "\n";

      if ((i % printingRate) == 0) cout << "\nGen = " << i;
   }

   cout << "\n gen = " << i-1;

   assign_rank_and_crowding_distance(parent_pop);
   report_params(file_params);
   report_pop(parent_pop, file_finalPop);
//   cerr << "\nParent pop size = " << parent_pop->getSize() << endl;
   report_feasible(parent_pop, file_feasiblePop);
   report_objSpace(parent_pop, file_objSpace);
   report_varSpace(parent_pop, file_varSpace);

   closeInfoFiles();
}


//string NSGA::rangesToStr(int digitsPrec, const char *delimiter) {
//   ostringstream strStream;
//
//   strStream << fixed;
//   for (int i = 0; i < nobj-1; ++i) {
//      strStream << setw(5 + digitsPrec) << setprecision(1 + digitsPrec)
//      << "[" << range_ndset.at(i).first << delimiter
//             << range_ndset.at(i).second << "]\n";
//   }
//
//   strStream << setw(5 + digitsPrec) << setprecision(1 + digitsPrec)
//             << "[" << range_ndset.at(nobj-1).first << delimiter
//                    << range_ndset.at(nobj-1).second << "]";
//
//   return strStream.str();
//}


void NSGA::openInfoFiles() {
   file_firstPop.open(name_firstPop.c_str(), ios::out);       //fpt1
   file_everyPop.open(name_everyPop.c_str(), ios::out);       //fpt4
   file_finalPop.open(name_finalPop.c_str(), ios::out);       //fpt2
   file_feasiblePop.open(name_feasiblePop.c_str(), ios::out); //fpt3
   file_params.open(name_params.c_str(), ios::out);           //ftp5

   file_everySubParentPop.open(name_everySubParentPop.c_str(), ios::out);
   file_everySubChildPop.open(name_everySubChildPop.c_str(), ios::out);
   file_varSpace.open(name_varSpace.c_str(), ios::out);
   file_objSpace.open(name_objSpace.c_str(), ios::out);
   file_inParentPop.open(name_inParentPop.c_str(), ios::out);
   file_subspaces.open(name_subspaces.c_str(), ios::out);
   file_subspacesInfo.open(name_subspacesInfo.c_str(), ios::out);
   file_onlineInd.open(name_onlineInd.c_str(), ios::out);

   file_intePhaseLength.open("intePhaseLenght.out", ios::out);
   file_partPhaseLength.open("partPhaseLenght.out", ios::out);
   file_convSubspaces.open("convSubspaces.out", ios::out);
   file_convWholeSpace.open("convWholeSpace.out", ios::out);
   file_parentShare.open("parentShare.out", ios::out);

   file_archive.open("ndominated_archive.out", ios::out);

   file_firstPop << "# This file contains the data of initial population.\n";
   file_finalPop << "# This file contains the data of final population.\n";
   //file_feasiblePop << "# This file contains the data of final feasible population (if found).\n";
   file_everyPop << "# This file contains the data of all generations.\n";

   ostringstream headerInfo(ostringstream::out);
   headerInfo << "# of objectives = " << nobj  << ", # of constraints = " << ncon
              << ", # of real_var = " << nreal << ", # of bits of bin_var = "
              << bitlength << "\n";

   file_firstPop << headerInfo.str();
   file_finalPop << headerInfo.str();
   file_everyPop << headerInfo.str();
   //file_feasiblePop << headerInfo.str();
}

void NSGA::closeInfoFiles() {
   file_firstPop.close();
   file_everyPop.close();
   file_finalPop.close();
   file_feasiblePop.close();
   file_params.close();

   file_everySubParentPop.close();
   file_everySubChildPop.close();
   file_varSpace.close();
   file_objSpace.close();
   file_inParentPop.close();
   file_subspaces.close();
   file_subspacesInfo.close();
   file_onlineInd.close();
   file_intePhaseLength.close();
   file_partPhaseLength.close();
   file_convSubspaces.close();
   file_convWholeSpace.close();
   file_parentShare.close();
   file_archive.close();
}


//void NSGA::subspaces_fill_nondominated_sort(
//         Population *mixed_pop,
//         Population *new_pop,
//         vector<vector<int> > &oSpaces)
//{
//   cout << "\n\n\nEntró al subspaces_fill\n" << endl;
//   /* Size of the subpopulation ranked in each objective subspace. */
//   int subSpacePopSize = parent_pop->getSize() / oSpaces.size(); /* integer division */
//
//   /* For each objective subspace */
//   int initPosition = 0; /* position in parent_pop to copy the whole sub_parent_pop */
//   for (unsigned s = 0; s < oSpaces.size(); ++s) {
//      /* In case of parent_pop's size was not divisible by numSubSpaces, the last
//         subpopulation will have the total pop - the individuals already distributed. */
//      if (s == (oSpaces.size() - 1)) {
////         cerr << "\n This is the last subpopulation..." << endl;
//         subSpacePopSize = parent_pop->getSize() - (s * subSpacePopSize);
//      }
//
//      Population subSpace_parent_pop(subSpacePopSize);
//      subSpace_parent_pop.cloneIndInitialize(mixed_pop->ind[0]);
//      //TODO: Inicializar la población.
////      cerr << "\nThe subpopulation " << s
////           << " has been initaliazed with " << subSpacePopSize << " individuals.";
//
//      /* Truncate mixed_pop by nondominated fronts. */
//      fill_nondominated_sort(mixed_pop, &subSpace_parent_pop, oSpaces[s]);
//      report_pop_objs(&subSpace_parent_pop, file_everySubParentPop);
//      file_everySubParentPop << "\n\n";
//      parent_pop->copyPopopulation(subSpace_parent_pop, initPosition);
//
//      //cout << "\nGlobal Parent pop\n";
//      //report_pop_objs(parent_pop, cout);
//
//      initPosition += subSpacePopSize;
//   }
//}

//void NSGA::subspaces_SortAndTruncationAndSelection(
//		Population *mixed_pop,
//		Population *new_parent_pop,
//		Population *new_child_pop,
//		vector<vector<int> > &oSpaces)
//{
//	report_pop_objs(mixed_pop, file_inParentPop);
//	file_inParentPop << "\n\n";
//
//   /* Size of the subpopulation ranked in each objective subspace. */
//   int subSpacePopSize = new_parent_pop->getSize() / oSpaces.size(); /* integer division */
//
//   /* For each objective subspace */
//   int initPosition = 0; /* position in child_pop to start the copy of sub_child_pop */
//   for (unsigned s = 0; s < oSpaces.size(); ++s) {
//      /* In case of parent_pop's size was not divisible by numSubSpaces, the last
//         subpopulation will have the total pop - the individuals already distributed. */
//      if (s == (oSpaces.size() - 1)) {
////         cerr << "\n This is the last subpopulation..." << endl;
//         subSpacePopSize = new_parent_pop->getSize() - (s * subSpacePopSize);
//      }
//
//      Population subSpace_parent_pop(subSpacePopSize);
//      Population subSpace_child_pop(subSpacePopSize);
//      subSpace_parent_pop.cloneIndInitialize(mixed_pop->ind[0]);
//      subSpace_child_pop.cloneIndInitialize(mixed_pop->ind[0]);
////      cerr << "\n\nThe child subpopulation has size " << subSpace_child_pop.size << " and is :\n";
////      report_pop_objs(&subSpace_child_pop, cout);
//
//      /* Truncate mixed_pop by nondominated fronts. */
//      fill_nondominated_sort(mixed_pop, &subSpace_parent_pop, oSpaces[s]);
//      //cout << "SubParentPop Size = " << subSpace_parent_pop.getSize() << endl;
//      report_pop_objs(&subSpace_parent_pop, file_everySubParentPop);
//      file_everySubParentPop << "\n\n";
//      new_parent_pop->copyPopopulation(subSpace_parent_pop, initPosition);
//
//      selection(&subSpace_parent_pop, &subSpace_child_pop);
//      //cout << "SubChildPop Size = " << subSpace_child_pop.getSize() << endl;
//      decode_pop(&subSpace_child_pop);
//      evaluate_pop(&subSpace_child_pop);
//      report_pop_objs(&subSpace_child_pop, file_everySubChildPop);
//      file_everySubChildPop << "\n\n";
//      new_child_pop->copyPopopulation(subSpace_child_pop, initPosition);
//
//      //cout << "\nGlobal Parent pop\n";
//      //report_pop_objs(parent_pop, cout);
//
//      initPosition += subSpacePopSize;
//   }
//}


/*
void NSGA::allocate_memory_pop (population *pop, int size) {
   int i;
   pop->ind = (individual *)malloc(size*sizeof(individual));
   for (i=0; i<size; i++)
   {
       allocate_memory_ind (&(pop->ind[i]));
   }
   return;
}
void NSGA::allocate_memory_ind (individual *ind) {
   int j;
   if (nreal != 0)
   {
       ind->xreal = (double *)malloc(nreal*sizeof(double));
   }
   if (nbin != 0)
   {
       ind->xbin = (double *)malloc(nbin*sizeof(double));
       ind->gene = (int **)malloc(nbin*sizeof(int));
       for (j=0; j<nbin; j++)
       {
           ind->gene[j] = (int *)malloc(nbits[j]*sizeof(int));
       }
   }
   ind->obj = (double *)malloc(nobj*sizeof(double));
   if (ncon != 0)
   {
       ind->constr = (double *)malloc(ncon*sizeof(double));
   }
   return;
}
void NSGA::deallocate_memory_pop (population *pop, int size) {
   int i;
   for (i=0; i<size; i++)
   {
       deallocate_memory_ind (&(pop->ind[i]));
   }
   free (pop->ind);
   return;
}
void NSGA::deallocate_memory_ind (individual *ind) {
   int j;
   if (nreal != 0)
   {
       free(ind->xreal);
   }
   if (nbin != 0)
   {
       for (j=0; j<nbin; j++)
       {
           free(ind->gene[j]);
       }
       free(ind->xbin);
       free(ind->gene);
   }
   free(ind->obj);
   if (ncon != 0)
   {
       free(ind->constr);
   }
   return;
}
*/

/*
double NSGA::maximum (double a, double b) {
   if ( a > b)
       return(a);

   return (b);
}

double NSGA::minimum (double a, double b) {
   if (a<b)
      return (a);

   return (b);
}
*/

void NSGA::crossover(Individual *parent1, Individual *parent2,
                     Individual *child1,  Individual *child2)
{
   if (nreal != 0)
       realcross(parent1, parent2, child1, child2);

   if (nbin != 0)
       bincross(parent1, parent2, child1, child2);

   return;
}

void NSGA::realcross(Individual *parent1, Individual *parent2,
                     Individual *child1,  Individual *child2)
{
   int i;
   double rand;
   double y1, y2, yl, yu;
   double c1, c2;
   double alpha, beta, betaq;

   if (r->randomperc() <= pcross_real)
   {
   	if (verbose) {
   	 cerr << "\nSe cruzarán" << endl;
   	}
       nrealcross++;
       for (i=0; i<nreal; i++)
       {
           if (r->randomperc() <= 0.5 )
           {
               if (fabs(parent1->xreal[i] - parent2->xreal[i]) > EPS)
               {
                   if (parent1->xreal[i] < parent2->xreal[i])
                   {
                       y1 = parent1->xreal[i];
                       y2 = parent2->xreal[i];
                   }
                   else
                   {
                       y1 = parent2->xreal[i];
                       y2 = parent1->xreal[i];
                   }
                   yl = range_realvar[i].first;
                   yu = range_realvar[i].second;
                   rand = r->randomperc();
                   beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                   alpha = 2.0 - pow(beta,-(eta_c+1.0));
                   if (rand <= (1.0/alpha))
                   {
                       betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                   }
                   else
                   {
                       betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                   }
                   c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                   beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                   alpha = 2.0 - pow(beta,-(eta_c+1.0));
                   if (rand <= (1.0/alpha))
                   {
                       betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                   }
                   else
                   {
                       betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                   }
                   c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                   if (c1<yl)
                       c1=yl;
                   if (c2<yl)
                       c2=yl;
                   if (c1>yu)
                       c1=yu;
                   if (c2>yu)
                       c2=yu;
                   //cerr << "\nSize xreal = " << child1->xreal.size() << endl;
                   //cerr << "\nvariable: " << i << endl;

                   if (r->randomperc() <= 0.5)
                   {
                       child1->xreal[i] = c2;
                       child2->xreal[i] = c1;
                   }
                   else
                   {
                       child1->xreal[i] = c1;
                       child2->xreal[i] = c2;
                   }
               }
               else
               {
                   child1->xreal[i] = parent1->xreal[i];
                   child2->xreal[i] = parent2->xreal[i];
               }
           }
           else
           {
               child1->xreal[i] = parent1->xreal[i];
               child2->xreal[i] = parent2->xreal[i];
           }
       }
   }
   else
   {
   	if (verbose) {
   	cerr << "\nSerán como los padres." << endl;
   	}
       for (i=0; i<nreal; i++)
       {
           child1->xreal[i] = parent1->xreal[i];
           child2->xreal[i] = parent2->xreal[i];
       }
   }

   if (verbose) {
      cout << "\npcx1 : ";
      parent1->printObjs(cout);
      cout << "\npcx2 : ";
      parent2->printObjs(cout);

   cout << "\ncx1 : ";
   child1->printObjs(cout);
   cout << "\ncx2 : ";
   child2->printObjs(cout);
   }
}

void NSGA::bincross(Individual *parent1, Individual *parent2, Individual *child1, Individual *child2)
{
   int i, j;
   double rand;
   int temp, site1, site2;
   for (i=0; i<nbin; i++)
   {
       rand = r->randomperc();
       if (rand <= pcross_bin)
       {
           nbincross++;
           site1 = r->rnd(0,nbits[i]-1);
           site2 = r->rnd(0,nbits[i]-1);
           if (site1 > site2)
           {
               temp = site1;
               site1 = site2;
               site2 = temp;
           }
           for (j=0; j<site1; j++)
           {
               child1->gene[i][j] = parent1->gene[i][j];
               child2->gene[i][j] = parent2->gene[i][j];
           }
           for (j=site1; j<site2; j++)
           {
               child1->gene[i][j] = parent2->gene[i][j];
               child2->gene[i][j] = parent1->gene[i][j];
           }
           for (j=site2; j<nbits[i]; j++)
           {
               child1->gene[i][j] = parent1->gene[i][j];
               child2->gene[i][j] = parent2->gene[i][j];
           }
       }
       else
       {
           for (j=0; j<nbits[i]; j++)
           {
               child1->gene[i][j] = parent1->gene[i][j];
               child2->gene[i][j] = parent2->gene[i][j];
           }
       }
   }
   return;
}

void NSGA::assign_crowding_distance_list(Population *pop, myList *lst, int front_size, vector<int> &objSet)
{
 //  int **obj_array;
 //  int *dist;
   unsigned i;
   int j;
   myList *temp;
   temp = lst;
   if (front_size==1)
   {
       pop->ind[lst->index]->crowd_dist = INF;
       return;
   }
   if (front_size==2)
   {
       pop->ind[lst->index]->crowd_dist = INF;
       pop->ind[lst->child->index]->crowd_dist = INF;
       return;
   }
   //obj_array = (int **) malloc(objSet.size() * sizeof(int *));
   vector<vector<int> > obj_array(objSet.size());
   vector<int> dist(front_size);
   //dist = (int *) malloc(front_size*sizeof(int));
   for (i=0; i< objSet.size(); i++)
       //obj_array[i] = (int *) malloc(front_size*sizeof(int));
      obj_array[i].resize(front_size);

   for (j=0; j<front_size; j++)
   {
       dist[j] = temp->index;
       temp = temp->child;
   }
   assign_crowding_distance(pop, dist, obj_array, objSet, front_size);
//   free (dist);
//   for (i=0; i<objSet.size(); i++)
//       free (obj_array[i]);

//   free (obj_array);
//   return;
}

void NSGA::assign_crowding_distance_list(Population *pop, myList *lst, int front_size)
{
   int **obj_array;
   int *dist;
   int i, j;
   myList *temp;
   temp = lst;
   if (front_size==1)
   {
       pop->ind[lst->index]->crowd_dist = INF;
       return;
   }
   if (front_size==2)
   {
       pop->ind[lst->index]->crowd_dist = INF;
       pop->ind[lst->child->index]->crowd_dist = INF;
       return;
   }
   obj_array = (int **) malloc(nobj * sizeof(int *));
   dist = (int *) malloc(front_size*sizeof(int));
   for (i=0; i< nobj; i++)
       obj_array[i] = (int *) malloc(front_size*sizeof(int));

   for (j=0; j<front_size; j++)
   {
       dist[j] = temp->index;
       temp = temp->child;
   }
   assign_crowding_distance(pop, dist, obj_array, front_size);
   free (dist);
   for (i=0; i<nobj; i++)
       free (obj_array[i]);

   free (obj_array);
   return;
}

void NSGA::assign_crowding_distance_indices(Population *pop, int c1, int c2, vector<int> &objSet)
{
//   int **obj_array;
//   int *dist;
   unsigned i;
   int j;
   int front_size;
   front_size = c2-c1+1;
   if (front_size==1)
   {
       pop->ind[c1]->crowd_dist = INF;
       return;
   }
   if (front_size==2)
   {
       pop->ind[c1]->crowd_dist = INF;
       pop->ind[c2]->crowd_dist = INF;
       return;
   }
   //obj_array = (int **) malloc(objSet.size() * sizeof(int *));
   vector<vector<int> > obj_array(objSet.size());
   vector<int> dist(front_size);
   //dist = (int *) malloc(front_size * sizeof(int));
   for (i=0; i < objSet.size(); i++)
      obj_array[i].resize(front_size);
            (int *) malloc(front_size*sizeof(int));
      //obj_array[i] = (int *) malloc(front_size*sizeof(int));

   for (j=0; j < front_size; j++)
       dist[j] = c1++;

   assign_crowding_distance(pop, dist, obj_array, objSet, front_size);

/*
   for (i=0; i < objSet.size(); i++)
       free (obj_array[i]);

   free (obj_array);
*/
   return;
}

void NSGA::assign_crowding_distance_indices(Population *pop, int c1, int c2)
{
   int **obj_array;
   int *dist;
   int i, j;
   int front_size;
   front_size = c2-c1+1;
   if (front_size==1)
   {
       pop->ind[c1]->crowd_dist = INF;
       return;
   }
   if (front_size==2)
   {
       pop->ind[c1]->crowd_dist = INF;
       pop->ind[c2]->crowd_dist = INF;
       return;
   }
   obj_array = (int **) malloc(nobj * sizeof(int *));
   dist = (int *) malloc(front_size * sizeof(int));
   for (i=0; i < nobj; i++)
       obj_array[i] = (int *) malloc(front_size*sizeof(int));

   for (j=0; j < front_size; j++)
       dist[j] = c1++;

   assign_crowding_distance(pop, dist, obj_array, front_size);
   free (dist);

   for (i=0; i < nobj; i++)
       free (obj_array[i]);

   free (obj_array);
   return;
}

void NSGA::assign_crowding_distance(Population *pop, vector<int> &dist,
                                    vector<vector<int> > &obj_array, vector<int> &objSet,
                                    int front_size)
{
   unsigned i;
   int j;

   for (i=0; i < objSet.size(); i++)
   {
       for (j=0; j < front_size; j++)
           obj_array[i][j] = dist[j];

       quicksort_front_obj(pop, objSet[i], obj_array[i], front_size);
   }

   for (j = 0; j < front_size; j++)
       pop->ind[ dist[j] ]->crowd_dist = 0.0;

   for (i=0; i < objSet.size(); i++)
       pop->ind[obj_array[ i ][0]]->crowd_dist = INF;

   for (i=0; i < objSet.size(); i++)
   {
       for (j=1; j<front_size-1; j++)
       {
           if (pop->ind[obj_array[ i ][j]]->crowd_dist != INF)
           {
               if (pop->ind[obj_array[i][front_size-1]]->obj[ objSet[i] ] ==
                   pop->ind[obj_array[i][0]]->obj[ objSet[i] ])
               {
                   pop->ind[obj_array[i][j]]->crowd_dist += 0.0;
               }
               else
               {
                   pop->ind[obj_array[i][j]]->crowd_dist +=
                         (pop->ind[obj_array[i][j+1]]->obj[ objSet[i] ] -
                          pop->ind[obj_array[i][j-1]]->obj[ objSet[i] ]) /
                         (pop->ind[obj_array[i][front_size-1]]->obj[ objSet[i] ] -
                          pop->ind[obj_array[i][0]]->obj[ objSet[i] ]);
               }
           }
       }
   }

   for (j=0; j<front_size; j++)
   {
       if (pop->ind[dist[j]]->crowd_dist != INF)
           pop->ind[dist[j]]->crowd_dist = (pop->ind[dist[j]]->crowd_dist) / objSet.size();
   }
}

void NSGA::assign_crowding_distance(Population *pop, int *dist, int **obj_array, int front_size)
{
   int i, j;
   for (i=0; i < nobj; i++)
   {
       for (j=0; j < front_size; j++)
           obj_array[ i ][j] = dist[j];

       quicksort_front_obj(pop, i, obj_array[i], front_size);
   }

   for (j = 0; j < front_size; j++)
       pop->ind[ dist[j] ]->crowd_dist = 0.0;

   for (i=0; i < nobj; i++)
       pop->ind[obj_array[ i ][0]]->crowd_dist = INF;

   for (i=0; i < nobj; i++)
   {
       for (j=1; j<front_size-1; j++)
       {
           if (pop->ind[obj_array[ i ][j]]->crowd_dist != INF)
           {
               if (pop->ind[obj_array[i][front_size-1]]->obj[ i ] ==
                   pop->ind[obj_array[i][0]]->obj[ i ])
               {
                   pop->ind[obj_array[i][j]]->crowd_dist += 0.0;
               }
               else
               {
                   pop->ind[obj_array[i][j]]->crowd_dist +=
                         (pop->ind[obj_array[i][j+1]]->obj[ i ] -
                          pop->ind[obj_array[i][j-1]]->obj[ i ]) /
                         (pop->ind[obj_array[i][front_size-1]]->obj[ i ] -
                          pop->ind[obj_array[i][0]]->obj[ i ]);
               }
           }
       }
   }

   for (j=0; j<front_size; j++)
   {
       if (pop->ind[dist[j]]->crowd_dist != INF)
           pop->ind[dist[j]]->crowd_dist = (pop->ind[dist[j]]->crowd_dist) / nobj;
   }
}

void NSGA::decode_pop(Population *pop) {
   int i;
   if (nbin!=0)
   {
       for (i=0; i<popsize; i++)
       {
          pop->ind[i]->decode();
          //decode_ind(pop->ind[i]);
       }
   }
}
/*
void NSGA::decode_ind(Individual *ind) {
   int j, k;
   int sum;
   if (nbin!=0)
   {
       for (j=0; j<nbin; j++)
       {
           sum=0;
           for (k=0; k<nbits[j]; k++)
           {
               if (ind->gene[j][k]==1)
               {
                   sum += pow(2,nbits[j]-1-k);
               }
           }
           ind->xbin[j] = range_binvar[j].first +
                          (double) sum*(range_binvar max_binvar[j] - min_binvar[j])/(double)(pow(2,nbits[j])-1);
       }
   }
   return;
}
*/

int NSGA::check_dominance(Individual *a, Individual *b, vector<int> &objSet) {
   int flag1;
   int flag2;
   flag1 = 0;
   flag2 = 0;

   if (a->constr_violation<0 && b->constr_violation<0)
   {
      if (a->constr_violation > b->constr_violation)
         return 1;
      else if (a->constr_violation < b->constr_violation)
         return -1;
      else
         return 0;
   }
   else if (a->constr_violation < 0 && b->constr_violation == 0)
      return -1;
   else
   {
      if (a->constr_violation == 0 && b->constr_violation <0)
         return 1;
      else
      {
         for (unsigned i=0; i < objSet.size(); i++)
         {
            if (a->obj[ objSet[i] ] < b->obj[ objSet[i] ])
               flag1 = 1;
            else if (a->obj[ objSet[i] ] > b->obj[ objSet[i] ])
               flag2 = 1;
         }

         if (flag1==1 && flag2==0)
            return 1;
         else if (flag1==0 && flag2==1)
            return -1;
         else
            return 0;
      }
   }
}

int NSGA::check_dominance (Individual *a, Individual *b) {
   int i;
   int flag1;
   int flag2;
   flag1 = 0;
   flag2 = 0;

   if (a->constr_violation<0 && b->constr_violation<0)
   {
      if (a->constr_violation > b->constr_violation)
         return 1;
      else if (a->constr_violation < b->constr_violation)
         return -1;
      else
         return 0;
   }
   else if (a->constr_violation < 0 && b->constr_violation == 0)
      return -1;
   else
   {
      if (a->constr_violation == 0 && b->constr_violation <0)
         return 1;
      else
      {
         for (i=0; i<nobj; i++)
         {
            if (a->obj[i] < b->obj[i])
               flag1 = 1;
            else if (a->obj[i] > b->obj[i])
               flag2 = 1;
         }

         if (flag1==1 && flag2==0)
            return 1;
         else if (flag1==0 && flag2==1)
            return -1;
         else
            return 0;
      }
   }
}

void NSGA::evaluate_pop(Population *pop) {
   pop->evaluate_pop(mop);
//   for (int i=0; i < popsize; ++i)
//      evaluate_ind(&(pop->ind[i]));
}

void NSGA::fill_nondominated_sort(Population *mixed_pop, Population *new_pop, vector<int> &objSet)
{
   int flag;
   int i, j;
   int end;
   int front_size;
   int archieve_size;
   int rank=1;
   myList *pool;
   myList *elite;
   myList *temp1, *temp2;
   pool = (myList *) malloc(sizeof(myList));
   elite = (myList *) malloc(sizeof(myList));
   front_size = 0;
   archieve_size=0;
   pool->index = -1;
   pool->parent = NULL;
   pool->child = NULL;
   elite->index = -1;
   elite->parent = NULL;
   elite->child = NULL;
   temp1 = pool;

//   cerr << "\n ******** entrada: fill_nondominated_sort ********" << endl;
   for (i=0; i< mixed_pop->getSize(); i++)
   {
       insert (temp1,i);
       temp1 = temp1->child;
   }
   i=0;
   do
   {
       temp1 = pool->child;
       insert (elite, temp1->index);
       front_size = 1;
       temp2 = elite->child;
       temp1 = del (temp1);
       temp1 = temp1->child;
       do
       {
           temp2 = elite->child;
           if (temp1==NULL)
           {
               break;
           }
           do
           {
               end = 0;
               //Revisar... Revisado.
               flag = check_dominance(mixed_pop->ind[temp1->index],
                                      mixed_pop->ind[temp2->index],
                                      objSet);
               if (flag == 1)
               {
                   insert(pool, temp2->index);
                   temp2 = del (temp2);
                   front_size--;
                   temp2 = temp2->child;
               }
               if (flag == 0)
                   temp2 = temp2->child;

               if (flag == -1)
                   end = 1;
           }
           while (end!=1 && temp2!=NULL);
           if (flag == 0 || flag == 1)
           {
               insert(elite, temp1->index);
               front_size++;
               temp1 = del(temp1);
           }
           temp1 = temp1->child;
       }
       while (temp1 != NULL);
       temp2 = elite->child;
       j=i;
/*
       cout << "arch_size = " << archieve_size << ", "
            << "fronto_size = " << front_size << ", "
            << "new_pop_size = " << new_pop->getSize() << endl;
*/
       if ( (archieve_size+front_size) <= new_pop->getSize())
       {
           do
           {
               //cout << "Escribiendo en posicion i = " << i;
               copy_ind(mixed_pop->ind[temp2->index], new_pop->ind[i]);
               //cout << "... éxito." << endl;
               new_pop->ind[i]->rank = rank;
               archieve_size += 1;
               temp2 = temp2->child;
               i+=1;
           } while (temp2 != NULL);

           //Revisar
           assign_crowding_distance_indices(new_pop, j, i-1, objSet);
           rank += 1;
       }
       else
       {
           //Revisar
           //cout << "\n Se truncará el último frente para que quepa.\n";
           crowding_fill(mixed_pop, new_pop, i, front_size, elite, objSet);
           archieve_size = new_pop->getSize();
           for (j=i; j < new_pop->getSize(); j++)
               new_pop->ind[j]->rank = rank;
       }
       temp2 = elite->child;
       do
       {
           temp2 = del (temp2);
           temp2 = temp2->child;
       }
       while (elite->child !=NULL);
   } while (archieve_size < new_pop->getSize());

   while (pool!=NULL)
   {
       temp1 = pool;
       pool = pool->child;
       free (temp1);
   }
   while (elite!=NULL)
   {
       temp1 = elite;
       elite = elite->child;
       free (temp1);
   }

//   cerr << "\n ******** salida: fill_nondominated_sort ********" << end;
}

void NSGA::fill_nondominated_sort(Population *mixed_pop, Population *new_pop)
{
   int flag;
   int i, j;
   int end;
   int front_size;
   int archieve_size;
   int rank=1;
   myList *pool;
   myList *elite;
   myList *temp1, *temp2;
   pool = (myList *) malloc(sizeof(myList));
   elite = (myList *) malloc(sizeof(myList));
   front_size = 0;
   archieve_size=0;
   pool->index = -1;
   pool->parent = NULL;
   pool->child = NULL;
   elite->index = -1;
   elite->parent = NULL;
   elite->child = NULL;
   temp1 = pool;

   for (i=0; i< mixed_pop->getSize(); i++)
   {
       insert (temp1,i);
       temp1 = temp1->child;
   }
   i=0;
   do
   {
       temp1 = pool->child;
       insert (elite, temp1->index);
       front_size = 1;
       temp2 = elite->child;
       temp1 = del (temp1);
       temp1 = temp1->child;
       do
       {
           temp2 = elite->child;
           if (temp1==NULL)
           {
               break;
           }
           do
           {
               end = 0;
               flag = check_dominance(mixed_pop->ind[temp1->index],
                                      mixed_pop->ind[temp2->index]);
               if (flag == 1)
               {
                   insert (pool, temp2->index);
                   temp2 = del (temp2);
                   front_size--;
                   temp2 = temp2->child;
               }
               if (flag == 0)
                   temp2 = temp2->child;

               if (flag == -1)
                   end = 1;
           }
           while (end!=1 && temp2!=NULL);
           if (flag == 0 || flag == 1)
           {
               insert(elite, temp1->index);
               front_size++;
               temp1 = del(temp1);
           }
           temp1 = temp1->child;
       }
       while (temp1 != NULL);
       temp2 = elite->child;
       j=i;
       if ( (archieve_size+front_size) <= new_pop->getSize())
       {
           do
           {
               copy_ind(mixed_pop->ind[temp2->index], new_pop->ind[i]);
               new_pop->ind[i]->rank = rank;
               archieve_size += 1;
               temp2 = temp2->child;
               i+=1;
           } while (temp2 != NULL);

           assign_crowding_distance_indices(new_pop, j, i-1);
           rank += 1;
       }
       else
       {
           crowding_fill(mixed_pop, new_pop, i, front_size, elite);
           archieve_size = new_pop->getSize();
           for (j=i; j < new_pop->getSize(); j++)
               new_pop->ind[j]->rank = rank;
       }
       temp2 = elite->child;
       do
       {
           temp2 = del (temp2);
           temp2 = temp2->child;
       }
       while (elite->child !=NULL);
   } while (archieve_size < new_pop->getSize());

   while (pool!=NULL)
   {
       temp1 = pool;
       pool = pool->child;
       free (temp1);
   }
   while (elite!=NULL)
   {
       temp1 = elite;
       elite = elite->child;
       free (temp1);
   }
}

void NSGA::crowding_fill(Population *mixed_pop, Population *new_pop, int count,
                         int front_size, myList *elite, vector<int> &objSet)
{
   myList *temp;
   int i, j;

   assign_crowding_distance_list(mixed_pop, elite->child, front_size, objSet);

   try{
      vector<int> dist(front_size);
      temp = elite->child;
      for (j=0; j < front_size; ++j)
      {
         dist.at(j) = temp->index;
         temp = temp->child;
      }
      quicksort_dist(mixed_pop, dist, front_size);

      for (i = count, j = front_size-1; i < new_pop->getSize(); ++i, --j)
         copy_ind(mixed_pop->ind[dist.at(j)], new_pop->ind[i]);

   }
   catch ( std::out_of_range outOfRange ) {// out_of_range exception
      cout << "\n\nException: " << outOfRange.what();
   }

}

void NSGA::crowding_fill(Population *mixed_pop, Population *new_pop, int count, int front_size, myList *elite)
{
   int *dist;
   myList *temp;
   int i, j;
   assign_crowding_distance_list(mixed_pop, elite->child, front_size);
   dist = (int *) malloc(front_size*sizeof(int));
   temp = elite->child;
   for (j=0; j<front_size; j++)
   {
       dist[j] = temp->index;
       temp = temp->child;
   }
   quicksort_dist(mixed_pop, dist, front_size);

   for (i=count, j=front_size-1; i<popsize; i++, j--)
       copy_ind(mixed_pop->ind[dist[j]], new_pop->ind[i]);

   free (dist);
}


void NSGA::insert(myList *node, int x) {
   myList *temp;
   if (node==NULL)
   {
       printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
       exit(1);
   }
   temp = (myList *) malloc(sizeof(myList));
   temp->index = x;
   temp->child = node->child;
   temp->parent = node;
   if (node->child != NULL)
   {
       node->child->parent = temp;
   }
   node->child = temp;
   return;
}

myList* NSGA::del(myList *node) {
   myList *temp;
   if (node==NULL)
   {
       printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
       exit(1);
   }
   temp = node->parent;
   temp->child = node->child;
   if (temp->child!=NULL)
   {
       temp->child->parent = temp;
   }
   free (node);
   return (temp);
}

void NSGA::merge(Population *pop1, Population *pop2, Population *pop3) {
   int i, k;
   for (i=0; i<popsize; i++)
       copy_ind(pop1->ind[i], pop3->ind[i]);

   for (i=0, k=popsize; i<popsize; i++, k++)
       copy_ind(pop2->ind[i], pop3->ind[k]);
}

void NSGA::copy_ind(Individual *ind1, Individual *ind2) {
   int i, j;
   ind2->rank = ind1->rank;
   ind2->constr_violation = ind1->constr_violation;
   ind2->crowd_dist = ind1->crowd_dist;
   if (nreal!=0)
   {
       for (i=0; i<nreal; i++)
           ind2->xreal[i] = ind1->xreal[i];
   }
   if (nbin!=0)
   {
       //cerr << "\nCopiando variables binarias...";
       for (i=0; i<nbin; i++)
       {
           ind2->xbin[i] = ind1->xbin[i];

           for (j=0; j<nbits[i]; j++)
               ind2->gene[i][j] = ind1->gene[i][j];
       }
   }

   //cerr << "\nCopiando objective values, copiando obj ";
   for (i=0; i<nobj; i++) {
       ind2->obj[i] = ind1->obj[i];
       //cerr << i << " ";
   }

   if (ncon!=0)
   {
    //  cerr << "\nCopiando valores de restricciones...";
       for (i=0; i<ncon; i++)
           ind2->constr[i] = ind1->constr[i];
   }
}

void NSGA::mutation_pop (Population *pop) {
   int i;

   for (i=0; i<popsize; i++)
       mutation_ind(pop->ind[i]);
}

void NSGA::mutation_ind(Individual *ind) {
   if (nreal != 0)
       real_mutate_ind(ind);

   if (nbin != 0)
       bin_mutate_ind(ind);
}

void NSGA::bin_mutate_ind (Individual *ind) {
   int j, k;
   double prob;
   for (j=0; j<nbin; j++)
   {
       for (k=0; k<nbits[j]; k++)
       {
           prob = r->randomperc();
           if (prob <=pmut_bin)
           {
               if (ind->gene[j][k] == 0)
                   ind->gene[j][k] = 1;
               else
                   ind->gene[j][k] = 0;

               nbinmut+=1;
           }
       }
   }
}

void NSGA::real_mutate_ind(Individual *ind) {
   int j;
   double rnd, delta1, delta2, mut_pow, deltaq;
   double y, yl, yu, val, xy;
   for (j=0; j<nreal; j++)
   {
       if (r->randomperc() <= pmut_real)
       {
           y = ind->xreal[j];
           yl = range_realvar[j].first;
           yu = range_realvar[j].second;
           delta1 = (y-yl)/(yu-yl);
           delta2 = (yu-y)/(yu-yl);
           rnd = r->randomperc();
           mut_pow = 1.0/(eta_m+1.0);
           if (rnd <= 0.5)
           {
               xy = 1.0-delta1;
               val = 2.0 * rnd+(1.0-2.0*rnd) * (pow(xy,(eta_m+1.0)));
               deltaq =  pow(val,mut_pow) - 1.0;
           }
           else
           {
               xy = 1.0-delta2;
               val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
               deltaq = 1.0 - (pow(val,mut_pow));
           }
           y = y + deltaq*(yu-yl);
           if (y<yl)
               y = yl;
           if (y>yu)
               y = yu;
           ind->xreal[j] = y;
           nrealmut += 1;
       }
   }
}

void NSGA::assign_rank_and_crowding_distance(Population *new_pop) {
   int flag;
   int i;
   int end;
   int front_size;
   int rank=1;
   myList *orig;
   myList *cur;
   myList *temp1, *temp2;
   orig = (myList *)malloc(sizeof(myList));
   cur = (myList *)malloc(sizeof(myList));
   front_size = 0;
   orig->index = -1;
   orig->parent = NULL;
   orig->child = NULL;
   cur->index = -1;
   cur->parent = NULL;
   cur->child = NULL;
   temp1 = orig;
   for (i=0; i<popsize; i++)
   {
       insert (temp1,i);
       temp1 = temp1->child;
   }
   do
   {
       if (orig->child->child == NULL)
       {
           new_pop->ind[orig->child->index]->rank = rank;
           new_pop->ind[orig->child->index]->crowd_dist = INF;
           break;
       }
       temp1 = orig->child;
       insert (cur, temp1->index);
       front_size = 1;
       temp2 = cur->child;
       temp1 = del (temp1);
       temp1 = temp1->child;
       do
       {
           temp2 = cur->child;
           do
           {
               end = 0;
               flag = check_dominance(new_pop->ind[temp1->index], new_pop->ind[temp2->index]);
               if (flag == 1)
               {
                   insert (orig, temp2->index);
                   temp2 = del (temp2);
                   front_size--;
                   temp2 = temp2->child;
               }
               if (flag == 0)
                   temp2 = temp2->child;

               if (flag == -1)
                   end = 1;

           } while (end != 1 && temp2 != NULL);

           if (flag == 0 || flag == 1)
           {
               insert(cur, temp1->index);
               front_size++;
               temp1 = del(temp1);
           }
           temp1 = temp1->child;

       } while (temp1 != NULL);

       temp2 = cur->child;

       do
       {
           new_pop->ind[temp2->index]->rank = rank;
           temp2 = temp2->child;

       } while (temp2 != NULL);

       assign_crowding_distance_list(new_pop, cur->child, front_size);

       temp2 = cur->child;
       do
       {
           temp2 = del(temp2);
           temp2 = temp2->child;
       }
       while (cur->child !=NULL);

       rank += 1;
   }
   while (orig->child != NULL);
   free (orig);
   free (cur);
}

void NSGA::quicksort_front_obj(Population *pop, int objcount, int obj_array[], int obj_array_size) {
   q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size-1);
}

void NSGA::quicksort_front_obj(Population *pop, int objcount,
                               vector<int> &obj_array, int obj_array_size)
{
   q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size-1);
}

void NSGA::q_sort_front_obj(Population *pop, int objcount,
                            vector<int> &obj_array, int left, int right)
{
   int index;
   int temp;
   int i, j;
   double pivot;
   if (left<right)
   {
       index = r->rnd(left, right);
       temp = obj_array[right];
       obj_array[right] = obj_array[index];
       obj_array[index] = temp;
       pivot = pop->ind[obj_array[right]]->obj[objcount];
       i = left-1;
       for (j=left; j<right; j++)
       {
           if (pop->ind[obj_array[j]]->obj[objcount] <= pivot)
           {
               i+=1;
               temp = obj_array[j];
               obj_array[j] = obj_array[i];
               obj_array[i] = temp;
           }
       }
       index=i+1;
       temp = obj_array[index];
       obj_array[index] = obj_array[right];
       obj_array[right] = temp;
       q_sort_front_obj (pop, objcount, obj_array, left, index-1);
       q_sort_front_obj (pop, objcount, obj_array, index+1, right);
   }
}

void NSGA::q_sort_front_obj(Population *pop, int objcount, int obj_array[], int left, int right) {
   int index;
   int temp;
   int i, j;
   double pivot;
   if (left<right)
   {
       index = r->rnd(left, right);
       temp = obj_array[right];
       obj_array[right] = obj_array[index];
       obj_array[index] = temp;
       pivot = pop->ind[obj_array[right]]->obj[objcount];
       i = left-1;
       for (j=left; j<right; j++)
       {
           if (pop->ind[obj_array[j]]->obj[objcount] <= pivot)
           {
               i+=1;
               temp = obj_array[j];
               obj_array[j] = obj_array[i];
               obj_array[i] = temp;
           }
       }
       index=i+1;
       temp = obj_array[index];
       obj_array[index] = obj_array[right];
       obj_array[right] = temp;
       q_sort_front_obj (pop, objcount, obj_array, left, index-1);
       q_sort_front_obj (pop, objcount, obj_array, index+1, right);
   }
}

void NSGA::quicksort_dist(Population *pop, vector<int> &dist, int front_size) {
   q_sort_dist(pop, dist, 0, front_size-1);
}

void NSGA::quicksort_dist(Population *pop, int *dist, int front_size) {
   q_sort_dist(pop, dist, 0, front_size-1);
}

void NSGA::q_sort_dist(Population *pop, vector<int> &dist, int left, int right) {
   int index;
   int temp;
   int i, j;
   double pivot;

   try {
      if (left<right)
      {
         index = r->rnd (left, right);
         temp = dist.at(right);
         dist.at(right) = dist.at(index);
         dist.at(index) = temp;
         pivot = pop->ind[dist.at(right)]->crowd_dist;
         i = left-1;
         for (j=left; j < right; j++)
         {
            if (pop->ind[dist.at(j)]->crowd_dist <= pivot)
            {
               i+=1;
               temp = dist.at(j);
               dist.at(j) = dist.at(i);
               dist.at(i) = temp;
            }
         }
         index=i+1;
         temp = dist.at(index);
         dist.at(index) = dist.at(right);
         dist.at(right) = temp;
         q_sort_dist(pop, dist, left, index-1);
         q_sort_dist(pop, dist, index+1, right);
      }
   } catch ( std::out_of_range outOfRange ) {// out_of_range exception
      cout << "\n\nException: " << outOfRange.what();
   }
}

void NSGA::q_sort_dist(Population *pop, int *dist, int left, int right) {
   int index;
   int temp;
   int i, j;
   double pivot;

   try {
      if (left<right)
      {
         index = r->rnd (left, right);
         temp = dist[right];
         dist[right] = dist[index];
         dist[index] = temp;
         pivot = pop->ind[dist[right]]->crowd_dist;
         i = left-1;
         for (j=left; j < right; j++)
         {
            if (pop->ind[dist[j]]->crowd_dist <= pivot)
            {
               i+=1;
               temp = dist[j];
               dist[j] = dist[i];
               dist[i] = temp;
            }
         }
         index=i+1;
         temp = dist[index];
         dist[index] = dist[right];
         dist[right] = temp;
         q_sort_dist(pop, dist, left, index-1);
         q_sort_dist(pop, dist, index+1, right);
      }
   } catch ( std::out_of_range outOfRange ) {// out_of_range exception
      cout << "\n\nException: " << outOfRange.what();
   }
}

void NSGA::selection(Population *old_pop, Population *new_pop) {
   int temp;
   int i;
   int rand;
   Individual *parent1, *parent2;
   vector<int> a1(old_pop->size);
   vector<int> a2(old_pop->size);
   //int *a1, *a2;
   //a1 = (int *) malloc(popsize*sizeof(int));
   //a2 = (int *) malloc(popsize*sizeof(int));

   for (i=0; i < old_pop->size; i++)
       a1[i] = a2[i] = i;

   for (i=0; i< old_pop->size; i++)
   {
       rand = r->rnd(i, old_pop->size-1);
       temp = a1[rand];
       a1[rand] = a1[i];
       a1[i] = temp;
       rand = r->rnd(i, old_pop->size-1);
       temp = a2[rand];
       a2[rand] = a2[i];
       a2[i] = temp;
   }

   for (i=0; i < old_pop->size; i+=4)
   {
   	 if (verbose) {
   		 cout << "\nTorneo " << i+1 << " entre " << endl;
          cout << "\np = " << a1[i] << " : ";
          old_pop->ind[a1[i]]->printObjs(cout);
          cout << "\np = " << a1[i+1] << " : ";
          old_pop->ind[a1[i+1]]->printObjs(cout);
          cout << "\ny ";
          cout << "\np = " << a1[i+2] << " : ";
          old_pop->ind[a1[i+2]]->printObjs(cout);
          cout << "\np = " << a1[i+3] << " : ";
          old_pop->ind[a1[i+3]]->printObjs(cout);
   	 }
       parent1 = tournament(old_pop->ind[a1[i]],   old_pop->ind[a1[i+1]]);
       parent2 = tournament(old_pop->ind[a1[i+2]], old_pop->ind[a1[i+3]]);
       if (verbose) {
          cout << "\n\nGanadores: " << "\n";
          parent1->printObjs(cout);
          cout << "\n";
          parent2->printObjs(cout);
       }
       crossover(parent1, parent2, new_pop->ind[i], new_pop->ind[i+1]);

       if (verbose) {
       cout << "\ni = " << i << " : ";
       new_pop->ind[i]->printObjs(cout);
       cout << "\ni = " << i+1 << " : ";
       new_pop->ind[i+1]->printObjs(cout);
       }

       if (verbose) {
   		 cout << "\nTorneo " << i+1 << " entre " << endl;
          cout << "\np = " << a2[i] << " : ";
          old_pop->ind[a2[i]]->printObjs(cout);
          cout << "\np = " << a2[i+1] << " : ";
          old_pop->ind[a2[i+1]]->printObjs(cout);
          cout << "\ny ";
          cout << "\np = " << a2[i+2] << " : ";
          old_pop->ind[a2[i+2]]->printObjs(cout);
          cout << "\np = " << a2[i+3] << " : ";
          old_pop->ind[a2[i+3]]->printObjs(cout);
       }

       parent1 = tournament(old_pop->ind[a2[i]], old_pop->ind[a2[i+1]]);
       parent2 = tournament(old_pop->ind[a2[i+2]], old_pop->ind[a2[i+3]]);
       if (verbose) {
          cout << "\n\nGanadores: " << "\n";
          parent1->printObjs(cout);
          cout << "\n";
          parent2->printObjs(cout);
       }

       crossover(parent1, parent2, new_pop->ind[i+2], new_pop->ind[i+3]);

       if (verbose) {
       cout << "\ni = " << i+2 << " : ";
       new_pop->ind[i+2]->printObjs(cout);
       cout << "\ni = " << i+3 << " : ";
       new_pop->ind[i+3]->printObjs(cout);
       }
   }
}

Individual* NSGA::tournament (Individual *ind1, Individual *ind2) {
   int flag;
   flag = check_dominance (ind1, ind2);
   if (flag==1)
   {
       return (ind1);
   }
   if (flag==-1)
   {
       return (ind2);
   }
   if (ind1->crowd_dist > ind2->crowd_dist)
   {
       return(ind1);
   }
   if (ind2->crowd_dist > ind1->crowd_dist)
   {
       return(ind2);
   }
   if ((r->randomperc()) <= 0.5)
   {
       return(ind1);
   }
   else
   {
       return(ind2);
   }
}

void NSGA::report_pop_objs(Population *pop, ostream &output) const
{
   int i, j;

   for (i=0; i < pop->getSize(); i++)
   {
       for (j=0; j < nobj; j++)
          output << fixed << setw(5 + 6) << setprecision(6)
                 << pop->ind[i]->obj[j];

       //output << setw(5) << pop->ind[i]->rank;

       //output << fixed << setw(5 + 4) << setprecision(4);

       //if (pop->ind[i]->crowd_dist >= INF)
       //   output << "INFTY";
       //else
       //   output << pop->ind[i]->crowd_dist;
       output << endl;
   }
}

/* Function to print the information of a population in a file */
void NSGA::report_pop(Population *pop, ostream &output) const
{
    int i, j, k;

    output << fixed << setw(5 + 8) << setprecision(8);
    for (i=0; i < pop->getSize(); i++)
    {
        for (j=0; j < nobj; j++)
           output << pop->ind[i]->obj[j] << "\t";

        if ( ncon!=0 )
        {
            for (j=0; j<ncon; j++)
               output << pop->ind[i]->constr[j] << "\t";
        }

        if ( nreal!=0 )
        {
            for (j=0; j<nreal; j++)
               output << pop->ind[i]->xreal[j] << "\t";
        }

        if ( nbin!=0 )
        {
            for (j=0; j<nbin; j++)
                for (k=0; k<nbits[j]; k++)
                   output << pop->ind[i]->gene[j][k] << "\t";
        }

        output << pop->ind[i]->constr_violation << "\t"
               << pop->ind[i]->rank << "\t"
               << pop->ind[i]->crowd_dist << endl;
    }
}

/* Function to print the information of feasible and non-dominated population in a file */
void NSGA::report_feasible(Population *pop, ostream &output) const
{
    int i, j, k;

    output << fixed << setw(5 + 8) << setprecision(8);

    for (i=0; i< pop->getSize(); i++)
    {
        if (pop->ind[i]->constr_violation == 0.0 && pop->ind[i]->rank == 1)
        {
            for (j=0; j<nobj; j++)
               output << pop->ind[i]->obj[j] << "\t";

            if (ncon != 0)  {
                for (j=0; j<ncon; j++)
                   output << pop->ind[i]->constr[j] << "\t";
            }
/*
            if (nreal != 0)  {
                for (j=0; j<nreal; j++)
                   output << pop->ind[i]->xreal[j] << "\t";
            }
*/
            if (nbin!=0)  {
                for (j=0; j<nbin; j++)
                    for (k=0; k<nbits[j]; k++)
                       output << pop->ind[i]->gene[j][k] << "\t";
            }
/*
            output << pop->ind[i]->constr_violation << "\t"
                   << pop->ind[i]->rank << "\t"
                   << pop->ind[i]->crowd_dist;
                   */

            output << endl;
        }
    }
}

/* Function to print the information of feasible and non-dominated population in a file */
void NSGA::report_objSpace(Population *pop, ostream &output) const
{
    int i, j;

    output << fixed << setw(5 + 8) << setprecision(8);

    for (i=0; i< pop->getSize(); i++)
    {
        if (pop->ind[i]->constr_violation == 0.0 && pop->ind[i]->rank == 1)
        {
            for (j=0; j<nobj-1; j++)
               output << pop->ind[i]->obj[j] << "\t";
            output << pop->ind[i]->obj[nobj-1];
            output << endl;
        }
    }
}

void NSGA::report_varSpace(Population *pop, ostream &output) const
{
    int i, j, k;

    output << fixed << setw(5 + 8) << setprecision(8);

    for (i=0; i< pop->getSize(); i++)
    {
        if (pop->ind[i]->constr_violation == 0.0 && pop->ind[i]->rank == 1)
        {
            if (nreal != 0)  {
                for (j=0; j<nreal; j++)
                   output << pop->ind[i]->xreal[j] << "\t";
            }

            if (nbin!=0)  {
                for (j=0; j<nbin; j++)
                    for (k=0; k<nbits[j]; k++)
                       output << pop->ind[i]->gene[j][k] << " ";
            }
            output << endl;
        }
    }
}

void NSGA::write_feasible(ostream &output) const
{
   report_feasible(parent_pop, output);
}

void NSGA::write_whole_pop(ostream &output) const
{
   report_pop(parent_pop, output);
}

void NSGA::report_params(ostream &output) const
{
   output << "# This file contains information about inputs as read by the program.\n";
   output << "\n Population size = " << popsize;
   output << "\n Number of generations = " << ngen;
   output << "\n Number of objective functions = " << nobj;
   output << "\n Number of constraints = " << ncon;
   output << "\n Number of real variables = " << nreal;

   if (nreal!=0)
   {
       for (int i=0; i< nreal; ++i)
       {
          output << "\n Lower limit of real variable " << i+1 << " = " << range_realvar[i].first;
          output << "\n Upper limit of real variable " << i+1 << " = " << range_realvar[i].second;
       }
       output << "\n Probability of crossover of real variable = " << pcross_real;
       output << "\n Probability of mutation of real variable = " << pmut_real;
       output << "\n Distribution index for crossover = " << eta_c;
       output << "\n Distribution index for mutation = " << eta_m;
   }
   output << "\n Number of binary variables = " << nbin;

   if (nbin!=0)
   {
       for (int i=0; i<nbin; ++i)
       {
          output << "\n Number of bits for binary variable " << i+1 << " = " << nbits[i];
          output << "\n Lower limit of binary variable " << i+1 << " = " << range_binvar[i].first;
          output << "\n Upper limit of binary variable " << i+1 << " = " << range_binvar[i].second;
       }
       output << "\n Probability of crossover of binary variable = " << pcross_bin;
       output << "\n Probability of mutation of binary variable = " << pmut_bin;
   }
   output << "\n Seed for random number generator = " << r->getSeed();

   if (nreal!=0)
   {
      output << "\n Number of crossover of real variable = " << nrealcross;
      output << "\n Number of mutation of real variable = " << nrealmut;
   }
   if (nbin!=0)
   {
      output << "\n Number of crossover of binary variable = " << nbincross;
      output << "\n Number of mutation of binary variable = " << nbinmut;
   }
}
