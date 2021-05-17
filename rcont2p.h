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


/* -----------------------------------------------------------------------
    Common declarations for rcont2 function.
----------------------------------------------------------------------- */

#include <assert.h>     /* Include standard assert function */
#include <limits.h>     /* Header file for implementation specific values */
#include <stdlib.h>     /* Header file for standard functions */
#include "fastexp2.h"   /* Header file for fastexp2 function */
#include "support.h"    /* Header file for support functions */


long double pangrand2(void);

/* -----------------------------------------------------------------------
    Restrictons:
    maxtot      Maximal total sum for the contingency table
----------------------------------------------------------------------- */

#define maxtot 100000   /* Maximal total sum */

/* -----------------------------------------------------------------------
    Debug pointers and memory allocation functions.

    PLEASE, READ THIS CAREFULLY BECAUSE THE FOLLOWING
    IS NOT COMMON PRACTICE FOR C/C++ PROGRAMMING!

    All rcont2 functions, except for rcont2_ access
    arrays thereby debug pointers instead of common
    C pointers.

    Debug pointer acts as common C pointer and also
    "knows" size of assigned array to perform bound
    check when arrays elements are accessed.
    If program attempt to access element outside
    array bounds, debug pointer calls prterr42(),

    Check of bounds requires C++ compiler with templates.
    To disable check of bounds either define macro NDEBUG
    or NCHECKPTR or use C compiler. If so, debug pointers
    will be implemented as common C pointers with no bound
    check (please, look for ANSI C implementation below).

    EXPTR(T) p;     - Declare debug pointer p to array
                      of type T.

    EXPTRTO(T,m,n)  - Construct debug pointer to existing
                      array m of n elements of type T.
                      Array may be represented thereby its
                      name, C pointer or debug pointer.
                      If the latter case, EXPTRTO macro
                      performs bound check.

    EXPTRCHECK(p,n) - Construct common C pointer for
                      debug pointer p and check if
                      elements 0..n-1 are accessible.

    EXPTRNEW(T,n)   - Allocate array for n elements of
                      type T and return its debug pointer.

    EXPTRDELETE(p)  - Deallocate array, allocated by
                      EXPTRNEW.
----------------------------------------------------------------------- */

/* #ifdef __cplusplus

/*
    Implementation for C++ with templates.
    Performs bound check if no NDEBUG or NCHECKPTR macro is defined
    Requires files exarray.h and exarray.c
*/

/* #include "exarray.h"    /* Header file for memory allocation functions */

/* #else /*__cplusplus*/

/*
    Implementation for ANSI C. No bound check.
    Uses custom sysalloc/sysfree functions from support package.
    Note, that EXPTRNEW allocated one extra item.
    This eliminates non-predictable consequences of common
    programming bug - indexing of non-existing element just
    after the last element. Bound check catches this error
*/

#define EXPTR(T)        T*
#define EXPTRTO(T,m,n)  (m)
#define EXPTRCHECK(p,n) (p)
#define EXPTRNEW(T,n)   ((T*) sysalloc ((n + 1), sizeof(T)))
#define EXPTRDELETE(p)  (sysfree (p))

/* #endif/*__cplusplus*/

/* -----------------------------------------------------------------------
    Constant values

    MAXTABLEROWS - Maximum number of rows in the contingency table.
    MAXTABLECOLS - Maximum number of cols in the contingency table.
    MAXTABLESUM  - Maximum sum of row/col in the contingency table.
----------------------------------------------------------------------- */

#define MAXTABLEROWS (SHRT_MAX / 2)  /* Suggest << SQRT (INT_MAX)       */
#define MAXTABLECOLS (SHRT_MAX / 2)  /* Suggest << SQRT (INT_MAX)       */
#define MAXTABLESUM  (SHRT_MAX / 2)  /* Suggest << SQRT (INT_MAX)       */

/* -----------------------------------------------------------------------
    Parameter passing conventions:

    1. Parameters type * represent variables, passed by
       rererence thereby C pointers.

       Input variables have const modifiers to prevent their
       accidental modification within function. Output and
       input/output variables have no const modifiers to allow
       their modification within function.

       For example: double *pre
                    const int nrow

       To pass a variable as a reference, use expression &name.

       For example:     &pre
                        &nrow

    2. Some functions accepts input variables as values.
       Appropriate parameter has no * sign and may have,
       or not to have const modifier.

       For example:     int size
                        int nitems

       To pass a value use name of variable or any expression:

       For example:     size
                        nrow + ncol

    3. Parameters markes as EXPTR(type) represent arrays,
       passed by their debug pointers.

       Input arrays have const modifiers to prevent their
       accidental modification within function. Output and
       input/output arrays have no const modifiers to allow
       their modification within function.

       For example:     EXPTR(double) fact
                        const EXPTR (int) irow

       To pass an array, use debug pointer of array or expression
       debug pointer + index.

       For example:     fact
                        irow + 1

       To produce debug pointer, allocate array thereby macro
       EXPTRNEW or apply macro EXPTRTO to existing array.

       For example:     EXPTR (double) fact = EXPTRNEW (double, 400)
                        static int rows [400];
                        EXPTR (int) irow = EXPTRTO (int, rows, 400)
                 
       Arrays, allocated by EXPTRNEW, are to be deallocated
       by means of EXPTRDELETE function.

       For example:     EXPTRDELETE (fact)
----------------------------------------------------------------------- */

/*
    Variables and functions defined in rcont2.c
*/

extern int   rcont2s (int nrow, int ncol,
                      int *nrowt,  int *ncolt,
                      int *matrix);


extern double rcont2f (void);

/*
    Function defined in rcont2_.c

    This function is converted by means of f2c
    and is included for the test purposes - its
    results are compared with result of manually
    converted rcount2 function.
*/

extern struct { float hop; } tempry_;

extern int  rcont2_ (int *nrow, int *ncol, int *nrowt, int *ncolt,
                     int *jwork, int *matrix, int *key, int *ifault);

/*
    Functions and variables defined in prterr.c
*/

extern void prcerr (int icode, const char *mes);
extern void prcerr1 (void);
extern void prcerr2 (void);
extern void prcerr3 (void);
extern void prcerr4 (void);
extern void prcerr5 (void);
extern void prcerr6 (void);

/* -----------------------------------------------------------------------
    Name:       RAND

    Purpose:    Generate random value in the range [0..1)

    Usage:      RAND

    Result:     The random value.
----------------------------------------------------------------------- */
#define RAND ((long double) pangrand2());

/*#define RAND ((float) rand() / (float) RAND_MAX)*/
