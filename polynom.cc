/************************************************************************
*                                                                       *
*               Structure Polynom and operations on it                  *
*                                                                       *
************************************************************************/

#include "includes.h"
#include "polynom.h"

extern "C" FILE *outfile;

namespace hyller {

/* init_Polynom(int n) -- initializes a general polynomial of nth power */

struct Polynom *init_Polynom(int n)
{
struct Polynom *polnm = NULL;
int i;

polnm = (struct Polynom*) malloc(sizeof(int) + sizeof(double*));
if (polnm == NULL)
 { printf("\tinit_Polynom : couldn't allocate memory. Bye\n");
   exit(0);
 }
else
 { (*polnm).power = n;
   (*polnm).body = (double*) init_array(n+1);
 }
return polnm;
}


/* destr_Polynom frees memory allocated for *P. Polynom and it's 
   body must have been allocated using init_Polynom */
void destr_Polynom(struct Polynom *P)
{ free((*P).body);
  free(P);
}

/******************** FUNCTIONS *********************************/

/* eval_Polynom gives a value of polynomial P in the point x */
double eval_Polynom(struct Polynom *P, double x)
{ 
int i;
double tmp=0;

if (x == 0)
   return (*P).body[0];
else
   { for(i = 0; i <= (*P).power; i++)
      tmp += (*P).body[i] * pow(x,i);
     return tmp;
   }
}

};

