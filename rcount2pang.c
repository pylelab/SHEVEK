/*******************************************************/
/******** 		Written By Phillip S. Pang      ********/
/*******************************************************/

/*******************************************************/
/******** MD/PhD Candidate, Columbia University ********/
/******** College of Physicians and Surgeons    ********/
/******** Dept. Of Biochemistry and Biophysics  ********/
/*******************************************************/

/*******************************************************/
/********    phillip.pang@stanfordalumni.org    ********/
/*******************************************************/


/*******************************************************/
/******** 		STATEMENT OF COPYRIGHT			********/

/*		    Copyright 2001 by The Trustees of          */
/*			Columbia University in the City of		   */
/*			New York. ALL RIGHTS RESERVED;			   */
			
/*******************************************************/

 

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*																	   */
/* THIS CODE WAS TRANSLATED, OPTIMIZED, and WRITTEN by RAUL Shakirov   */
/* from FORTRAN TO C, As a work-for-hire by request of Phillip S. Pang */
/*																	   */
/* Contact information for Raul: http://www.imach.uran.ru/rns/		   */
/*																	   */
/*																	   */
/* It has been further modified by Phillip S. Pang					   */
/*																	   */
/***********************************************************************/
/***********************************************************************/
/***********************************************************************/



/****************************************************************/
/****************************************************************/
/****************************************************************/
/*																*/
/*	The algorithms found within this file may be derivatives	*/
/*  of source code obtained from:								*/
/*    ALGORITHM AS 159  APPL. STATIST. (1981) VOL.30, NO.1		*/
/*	See: Patefield, Applied Statistics 30:91-97 (1981) 			*/
/*																*/
/****************************************************************/
/****************************************************************/
/****************************************************************/




#include <float.h>      /* Header file for standard math functions */
#include <math.h>       /* Header file for standard math functions */
#include "rcont2p.h"     /* Header file for rcont2 functions */

/* -----------------------------------------------------------------------
    ALGORITHM AS 159  APPL. STATIST. (1981) VOL.30, NO.1

    Name:       rcont2s
                rcont2f

    Purpose:    Generate random two-way table with given marginal totals

    Usage:      srand ((unsigned) time (NULL))
                rcont2s (nrow, ncol, nrowt, ncolt, matirix)                
                rcont2f()

       nrow   - The number of rows in the contingecy table.       (Input)
       ncol   - The number of columns in the contingency table.   (Input)
       nrowt  - Vector of length nrow contining row totals
                of the observed table                             (Input)
       ncolt  - Vector of length ncol contining column totals
                of the observed table                             (Input)
       matrix - nrow by ncol matrix for the contingency table
                to be generated                                  (Output)

    Notes:      srand   Standard C function which is used to set
                        pseudo-random seed (requires stdlib.h).

                time    Standard C function which is used to get
                        current time (requires time.h).

                rcont2s stores arguments and returns error code


                rcont2f does same as rcont2 but uses fast
                        table-driven exp function with 20 bit
                        presizion instead of standard exp function.
----------------------------------------------------------------------- */




/**************************************************************************/
/**************************************************************************/

/*
    Static variables and arrays are kept between calls of
    function because they are allocated and inilialized once
    when the program is loaded. Default initial value - 0.
*/

/*  Static variables for arguments. */

static int nrowm;               /* The number of rows - 1 */
static int ncolm;               /* The number of cols - 1 */
static int *nrowt;        /* Vector of row totals */
static int *ncolt;        /* Vector of col totals */
static int *matrix;       /* Matrix for the contingency table */
static int *jwork;        /* Last row of matrix */
static int ntotal;              /* Total sum */

/*
    Static work array for optimization purposes
    log - log-factorials (log(1) + log(2) + ... log(n))
    div - dividors       (1/n)
*/

static int nfill = 0;           /* Index of filled part of fact */
                                /* Work array */
static struct {double log, div; } fact [maxtot + 1];


/**************************************************************************/
/**************************************************************************/






/*************************************************************************/
/* -----------------------------------------------------------------------
    Name:       rcont2s

    Purpose:    Prepare to construct random matrix.

    Usage:      Ref. top of file.

    Returns:    Error code:
                   1 Number of rows must be greater than 1.
                   2 Number of columns must be greater than 1.
                   3 Row totals must be greater than 0.
                   4 Column totals must be greater than 0.
                   5 Sample size exceeds 5000.
                   6 Sum of row and column totals must be the same

    External referencies:
       log, prcerr1..6
----------------------------------------------------------------------- */
/************************************************************************/




int rcont2s (int nrow, int ncol,
             int *_nrowt,  int *_ncolt,
             int *_matrix)
{
    /* Local variables */

    int ntot = 0;

    /* Check for faults and prepare for future calls of rcont2x */

    if ((nrowm = nrow - 1) <= 0)
    {
        return 1;                       /* Error */
    }

    if ((ncolm = ncol - 1) <= 0)
    {
        return 2;                       /* Error */
    }

    /* Store arrays. */

    nrowt  = _nrowt;
    ncolt  = _ncolt;
    matrix = _matrix;

    /* EXPTRTO macro calculates debug pointer to sub-array. */

    jwork  = EXPTRTO (int, matrix + nrowm * ncol, ncol);

    /* Compute total sum */

    {
        int i, ilim;                /* Loop index and limit */
        ilim = nrow;
        for (i = 0; i < ilim; ++i)
        {
            int n = nrowt [i];
            if (n <= 0)
            {
                return 3;               /* Error */
            }
            ntot += n;
        }
    }

    if (ntot > maxtot)
    {
       return 5;                        /* Error */
    }

    /* Calculate log-factorials and dividors */

    {
        /* Local variables for this blick */

        int n = nfill;                  /* Loop index */
        double f = fact [n].log;         /* Fact value */

        for (++n; n <= ntot; ++n)
        {
            f += (double) log ((double) n);
            fact [n].log = f;
            fact [n].div = (double) 1./((double) n);
        }

        nfill = ntot;
    }

    /* Store total sum */

    ntotal = ntot;

    /* Check if row and column sums are the same */

    {
        int j, jlim;                /* Loop index and limit */
        jlim = ncol;
        for (j = 0; j < jlim; ++j)
        {
            int n = ncolt [j];
            if (n <= 0)
            {
                return 4;               /* Error */
            }
            ntot -= n;
        }
    }

    if (ntot != 0)
    {
       return 6;                        /* Error */
    }

    return 0;                           /* No errors */
}




/************************************************************************/
/* -----------------------------------------------------------------------
    Name:       rcont2f

    Purpose:    Construct random matrix using table-driven
                20 bit exp function.

    Usage:      Ref. top of file.

    Returns:    Probability of observing the sample table
                as double value.

    External referencies:
       rand, fastexp2
----------------------------------------------------------------------- */
/************************************************************************/



double rcont2f (void)

#define FASTEXP(value) fastexp2(value)  /* Use table-driven exp function */


{
    /* Local variables */

    int i, j;                   /* Loop index and limit */
    EXPTR(int) matp = matrix;   /* Current element to be calculated */

    double hop = 1.;             /* Probabilty value */

    int jc = ntotal;            /* Total sum for subsequent rows */
    int ib;                     /* Remainder of total sum */

    /* Copy ncolt to last row of matrix (jwork) */

    for (j = 0; j < ncolm; ++j)
    {
       jwork [j] = ncolt [j];
    }

    /* Cycle for rows */

    for (i = 0; i < nrowm; ++i)
    {
        /* Local variables for this cycle */

        int ic = jc;            /* Remainder of total sum */
        int ia = nrowt [i];     /* Remainder of row total */

        /* Update jc */

        jc -= ia;

        /* Cycle for elements of the row */

        for (j = 0; j < ncolm; ++j)
        {
            /*  Local variable for this cycle */

            int nlm;            /* Value of element to be computed */

            /* Test for zero entries in matrix */

            if (ic != 0)
            {
                /* Generate pseudo-random number */

                long double dummy; 

                /* Update remainders */

                const int ie = ic;          /* Remainder of total sum */
                const int id = jwork [j];   /* Remainder of column total */
                const int ii = (ib = ic - ia) - id;

                /* Calulate initial value of nlm */

                double fnlm = (double) ia * (double) id * fact [ie].div 
                        #if defined(_MSC_VER) && defined(_M_IX86)
                           * (double) (1. + FLT_EPSILON);
                        #else
                           + (double) 0.5;
                        #endif

                /* Update ic */

                ic -= id;

				dummy = (long double) RAND;

                /* Compute conditional expected value of matrix [i] [j] */

                for (;;)
                {
                    /* Local variables for this cycle */

                    double x, y, sumprb;       /* Work values */
                    int nll;                  /* nlm, nlm-1, ... */

                #if defined(_MSC_VER) && defined(_M_IX86)
                    __asm fld   fnlm
                    __asm fistp nlm
                #else 
                    nlm = (int) fnlm;
                #endif

                    x = (double) FASTEXP (
                            fact [ia].log + fact [ib].log +
                            fact [ic].log + fact [id].log - fact [ie].log -
                            fact [id - nlm].log - fact [ia - nlm].log -
                            fact [ii + nlm].log - fact [nlm].log
                            );

                    if (x >= dummy)
                    {
                        hop *= x;
                        break /* computing for */;
                    }

                    y = x;
                    sumprb = x;
                    nll = nlm;

                    for (;;)
                    {
                        /* Local variables for this cycle */

                        int lsp, lsm;           /* Work values */

                        if ((lsp = id - nlm) != 0 &&
                            (lsm = ia - nlm) != 0)
                        {
                            /* Increment entry in row i, column j */

                            ++nlm;

                            x *= ((double) lsp * fact [nlm].div) *
                                 ((double) lsm * fact [ii + nlm].div);

                            if ((sumprb += x) >= dummy)
                            {
                                hop *= x;
                                goto MATP /* break computing for */;
                            }
                        }
                        else
                        {
                            while ((lsp = nll) != 0 &&
                                   (lsm = ii + nll) != 0)
                            {
                                /* Decrement entry in row i, column j */

                                --nll;

                                y *= ((double) lsp * fact [id - nll].div) *
                                     ((double) lsm * fact [ia - nll].div);

                                if ((sumprb += y) >= dummy)
                                {
                                    hop *= y;
                                    nlm = nll;
                                    goto MATP /* break computing for */;
                                }
                            }

                            break /* for */;
                        }

                        if ((lsp = nll) != 0 &&
                            (lsm = ii + nll) != 0)
                        {
                            /* Decrement entry in row i, column j */

                            --nll;

                            y *= ((double) lsp * fact [id - nll].div) *
                                 ((double) lsm * fact [ia - nll].div);

                            if ((sumprb += y) >= dummy)
                            {
                                hop *= y;
                                nlm = nll;
                                goto MATP /* break computing for */;
                            }
                        }
                        else
                        {
                            while ((lsp = id - nlm) != 0 &&
                                   (lsm = ia - nlm) != 0)
                            {
                                /* Increment entry in row i, column j */

                                ++nlm;

                                x *= ((double) lsp * fact [nlm].div) *
                                     ((double) lsm * fact [ii + nlm].div);

                                if ((sumprb += x) >= dummy)
                                {                              
                                    hop *= x;
                                    goto MATP /* break computing for */;
                                }
                            }
                            break /* for */;
                        }
                    }

                    /* Generate new pseudo-random number */

                    dummy = (long double) sumprb*RAND;
                }

            MATP:
                /* Set element of matrix */

                *matp++ = nlm;
                jwork [j] -= nlm;
                ia -= nlm;
            }
            else /* ic == 0 */
            {
                /* All remaining elements in the row are to be zero */

                for (; j < ncolm; ++j)
                {
                    *matp++ = 0;
                }

                /*
                   Probably, it is sufficient to set ib = 0,
                   because ia should be 0
                */

                ib = ic - ia;
                ia = 0;

                break /* for */;
            }
        }

        /* Set last element in the row */

        *matp++ = ia;
    }

    /* Compute last element in last row of matrix */

    jwork [j] = ib - jwork [j - 1];

    return hop;
}



#undef  FASTEXP                         /* Undef fastexp */
