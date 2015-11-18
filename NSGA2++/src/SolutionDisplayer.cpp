/*
 * SolutionDisplayer.cpp
 *
 *  Created on: 11/11/2009
 *      Author: antonio
 */

#include <iostream>
#include "SolutionDisplayer.h"

using namespace std;

/* Function to print the information of a population in a file */
void SolutionDisplayer::report_pop(Population *pop, ostream &output)
{
   int i, j, k;

   for (i=0; i < pop->size; i++)
   {
       for (j=0; j < pop->nobj; j++)
          output << pop->ind[i]->obj[j] << "\t";

       if ( pop->ncon!=0 )
       {
           for (j=0; j < pop->ncon; j++)
              output << pop->ind[i]->constr[j] << "\t";
       }

       if ( pop->nreal!=0 )
       {
           for (j=0; j< pop->nreal; j++)
              output << pop->ind[i]->xreal[j] << "\t";
       }

       if ( pop->nbin!=0 )
       {
           for (j=0; j < pop->nbin; j++)
               for (k=0; k< pop->nbits[j]; k++)
                  output << pop->ind[i]->gene[j][k] << "\t";
       }

       cout << pop->ind[i]->constr_violation << "\t"
            << pop->ind[i]->rank << "\t"
            << pop->ind[i]->crowd_dist << endl;
   }
}



/* Function to print the information of feasible and non-dominated population in a file */
void SolutionDisplayer::report_feasible(Population *pop, ostream &output)
{
   int i, j, k;

   for (i=0; i< pop->size; i++)
   {
       if (pop->ind[i]->constr_violation == 0.0 && pop->ind[i]->rank == 1)
       {
           for (j=0; j< pop->nobj; j++)
              output << pop->ind[i]->obj[j] << "\t";

           if (pop->ncon != 0)  {
               for (j=0; j< pop->ncon; j++)
                  output << pop->ind[i]->constr[j] << "\t";
           }

           if (pop->nreal != 0)  {
               for (j=0; j< pop->nreal; j++)
                  output << pop->ind[i]->xreal[j] << "\t";
           }

           if (pop->nbin!=0)  {
               for (j=0; j<pop->nbin; j++)
                   for (k=0; k<pop->nbits[j]; k++)
                      output << pop->ind[i]->gene[j][k] << "\t";
           }

           output << pop->ind[i]->constr_violation << "\t"
                  << pop->ind[i]->rank << "\t"
                  << pop->ind[i]->crowd_dist << endl;
       }
   }
}

