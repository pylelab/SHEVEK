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



#include <math.h>       /* Header file for standard math functions */
#include <stdio.h>      /* Header file for standard io functions */
#include "fastexp2.h"   /* Header file for fastexp2 function */

/* -----------------------------------------------------------------------
    Key constants (ref. explanations for function fastexp2 below).
----------------------------------------------------------------------- */
#define EXPMINVAL (0)
#define EXPMAXVAL (8)
#define EXPBITS1  (9)
#define EXPBITS2  (11)

/* -----------------------------------------------------------------------
    Derivative constants 
    MULT1      pow (2, EXPBITS1)
    MULT2      pow (2, EXPBITS2)
    EXPTABLE1  num of elements in exptable1 (ref. fastexp2 below)
    EXPTABLE2  num of elements in exptable2 (ref. fastexp2 below)
    nmult12    multiplicator for fastexp2 
----------------------------------------------------------------------- */
#define MULT1     (2 << (EXPBITS1-1))
#define MULT2     (2 << (EXPBITS2-1))
#define EXPTABLE1 (MULT1 * EXPMAXVAL)
#define EXPTABLE2 (MULT2)
const float nmult12 = (- MULT1 * MULT2);

/* -----------------------------------------------------------------------
    Include file with content of tables exptable1 and exptable2.
    This file is built by call of fastexp2inc (ref. below).
----------------------------------------------------------------------- */
#include "fastexp2.inc"

/* -----------------------------------------------------------------------
    Name:       fastexp2ini

    Purpose:    Build tables for fastexp2 for given key constants.

    Usage:      fastexp2ini()

    Note:       Used for development purposes only!
----------------------------------------------------------------------- */

void fastexp2ini (void)
{
    int i;

    for (i = 0; i < EXPTABLE1; i++)
    {
        exptable1 [i] = (float)(exp (-((double)i) / MULT1));
    }

    for (i = 0; i < EXPTABLE2; i++)
    {
#if defined(_MSC_VER) && defined(_M_IX86)

        /*
            fastexp2 for VC++ and x86 processors uses coversion of
            double to nearest int (round mode) so we haven't add 0.5
            to exp argument here.
        */
        exptable2 [i] = (float)(exp (-((double)i) / (MULT1 * MULT2)));

#else /*_MSC_VER && _M_IX86*/

        /*
            fastexp2 for generic processors uses conversion of double
            to lower nearest int (floor mode) so we have add 0.5
            to exp argument for better approximation.
        */
        exptable2 [i] = (float)(exp (-((double)i+0.5) / (MULT1 * MULT2)));

#endif/*_MSC_VER && _M_IX86*/
    }
}

/* -----------------------------------------------------------------------
    Name:       fastexp2inc

    Purpose:    Print tables for fastexp2 for given key contatnts
                (file fastexp2.inc).

    Usage:      fastexp2inc()

    Note:       Used for development purposes only!
----------------------------------------------------------------------- */

void fastexp2inc (void)
{
    int i;

    printf ("/*\n"
        "    File with fastexp2 tables for case\n"
        "    EXPTABLE1 = %d, MULT1 = %d, EXPTABLE2 = %d, MULT2 = %d\n"
        "    The file was generated by function fastexp2inc.\n*/\n",
        EXPTABLE1, MULT1, EXPTABLE2, MULT2);

    printf ("\n#if EXPTABLE1 != %d || MULT1 != %d || "
                  "EXPTABLE2 != %d || MULT2 != %d\n",
                   EXPTABLE1, MULT1, EXPTABLE2, MULT2);

    printf ("\nstatic float exptable1 [EXPTABLE1];");

    printf ("\nstatic float exptable2 [EXPTABLE2];");

    printf ("\n\n#else /*EXPTABLE && MULT*/\n");

    printf ("\nstatic float exptable1 [EXPTABLE1] =\n{");

    for (i = 0;;)
    {
        if (i % 4 == 0) printf ("\n   ");
        printf (" %.7ef",
                (float)(exp (-((double)i) / MULT1))
               );
        if (++i < EXPTABLE1) printf (","); else break;
    }

    printf ("\n};\n");

    printf ("\nstatic float exptable2 [EXPTABLE2] =\n{");

    printf ("\n#if defined(_MSC_VER) && defined(_M_IX86)\n");

    for (i = 0;;)
    {
        if (i % 4 == 0) printf ("\n   ");
        /*
            fastexp2 for VC++ and x86 processors uses coversion of
            double to nearest int (round mode) so we haven't add 0.5
            to exp argument here.
        */
        printf (" %.7ef",
                (float)(exp (-((double)i) / (MULT1 * MULT2)))
               );
        if (++i < EXPTABLE2) printf (","); else break;
    }

    printf ("\n\n#else /*_MSC_VER && _M_IX86*/\n");

    for (i = 0;;)
    {
        if (i % 4 == 0) printf ("\n   ");
        /*
            fastexp2 for generic processors uses conversion of double
            to lower nearest int (floor mode) so we have add 0.5
            to exp argument for better approximation.
        */
        printf (" %.7ef",
                (float)(exp (-((double)i+0.5) / (MULT1 * MULT2)))
               );
        if (++i < EXPTABLE2) printf (","); else break;
    }

    printf ("\n\n#endif/*_MSC_VER && _M_IX86*/\n};\n");

    printf ("\n#endif/*EXPTABLE && MULT*/\n");
}

/* -----------------------------------------------------------------------
    Name:       fastexp2

    Purpose:    Fast table-driven exponent algorithm.

    Usage:      fastexp2 (value)

    Result:     The exp of value.

    Method:     exp (i1 + i2) = exp (i1) * exp (i2)

                i1 is integer part of value and first EXPBITS1 bits
                of fractional part and i2 is successive EXPBITS2 bits
                of fractional part.

    Function gets exp (i1) from array exptable1 and exp (i2) from
    array exptable2 for values between -EXPMAXVAL and -EXPMINVAL.
    Max relative error is 1 / pow (2, EXPBITS1 + EXPBITS2 + 1).
    Function returns EXPBITS1 + EXPBITS2 + 1 significant bits.

    If value > -EXMINVAL or <= -EXMAXVAL, standard exp function is called.
----------------------------------------------------------------------- */

double fastexp2 (double value)
{
    register int i1, i2;                /* Components of value */

    /* Extract valid bits to integer variable */

#if defined(_MSC_VER) && defined(_M_IX86)

    /*
        fastexp2 for VC++ and x86 processors uses coversion of
        double to nearest int (round mode).
    */

    __asm fld   nmult12
    __asm fmul  value
    __asm fistp i1

#else /*_MSC_VER && _M_IX86*/

    /*
        fastexp2 for generic processors uses conversion of double
        to lower nearest int (floor mode).
        On x86 it's less efficient then conversion to nearest int.        
    */

    i1 = (int) (value * nmult12);

#endif/*_MSC_VER && _M_IX86*/

    /* Divide valid bits to high and low parts */

    i2 = i1 & EXPTABLE2 - 1;            /* Get low  BITS2 bits */
    i1 >>= EXPBITS2;                    /* Get high BITS1 bits */

    if (i1 >= MULT1 * EXPMINVAL)
    {
        if (i1 < MULT1 * EXPMAXVAL)     /* EXPTABLE1 */
        {
            /* Use tables to get exp (i1) and exp (i2) */

            return (exptable1 [i1] * exptable2 [i2]);
        }
    }

    /* Use standard exp function. */

    return (exp (value));
}
