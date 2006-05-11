
#ifndef _hyller_polynom_h_
#define _hyller_polynom_h_

namespace hyller {

/************************************************************************
*									*
*		Structure Polynom and operations on it			*
*									*
************************************************************************/
/* Detailed descriptions are in the source file (polynom.c) */

/* Polynomial is a structure */
struct Polynom {
        int power;
        double *body;
        };
                        
/* Various initialization routines */
struct Polynom *init_Polynom(int size);
void destr_Polynom(struct Polynom *P);


/* Various functions on polynomials */
double eval_Polynom(struct Polynom *P, double x);

};

#endif
