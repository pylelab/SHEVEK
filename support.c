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



#include <malloc.h>     /* Header file for standard memory functions */
#include <stdio.h>      /* Header file for standard io functions */
#include <stdlib.h>     /* Header file for standard functions */
#include "support.h"    /* Header file for support functions */

/* -----------------------------------------------------------------------
   Static variable to check if memory heap is properly deallocated.
----------------------------------------------------------------------- */

static int syscount;

/* -----------------------------------------------------------------------
    Name:       sysalloc

    Purpose:    Allocate memory for array.

    Usage:      sysalloc (nitems, size)

    Arguments:
       nitems - Number of intes in array
       size   - size of item

    Result:   - Pointer to allocated array

    External referencies:
       syserr40
----------------------------------------------------------------------- */
 
void* sysalloc (int nitems, int size)
{
   void *p = calloc (nitems, size);

   if (p == NULL) syserr40();

   ++syscount;

   return p;
}

/* -----------------------------------------------------------------------
    Name:       sysfree

    Purpose:    Deallocate memory, allocated by sysalloc,

    Usage:      sysfree (p)

    Arguments:
       p      - Pointer to deallocated array
----------------------------------------------------------------------- */

void sysfree (void *p)
{
   if (p != NULL)
   {
      --syscount;

      free (p);
   }
}

/* -----------------------------------------------------------------------
    Name:       syschk

    Purpose:    Check if memory is properly deallocated.

    Usage:      syschk()
----------------------------------------------------------------------- */

void syschk (void)
{
   if (syscount)                        /* Print diagnostic to console */
   {
      fprintf (stderr, "\nSYSTEM WARNING: Memory leak detected.\n"
                       "It's an internal error. Ask developers.\n");
   }
}

/* -----------------------------------------------------------------------
    Name:       syserr

    Purpose:    Print an error message to console and exit.

    Usage:      prterr (icode, mes)

    Arguments:
      icode   - Int code for the error message.                 (Input)
      mes     - Character string containing the error message.  (Input)
----------------------------------------------------------------------- */

void syserr (int icode, const char *mes)
{
    fprintf (stderr, "\nSYSTEM ERROR: %d %s\n", icode, mes);
                                /* Print diagnostic to console */
    exit (2);                   /* Close all files and exit with code 1 */
}

/* -----------------------------------------------------------------------
    Name:       syserrN

    Purpose:    Print error message N to console and exit.

    Usage:      syserrN()
----------------------------------------------------------------------- */

void syserr40 (void)
{
    syserr (40, "Out of workspace.\n"
                "Select system with more RAM available");
}

void syserr41 (void)
{
    syserr (41, "Automatic expansion of array is not permissible.\n"
                "It's an internal error. Ask developers.");
}

void syserr42 (void)
{
    syserr (42, "Index is out of range.\n" 
                "It's an internal error. Ask developers.");
}

void syserr43 (void)
{
    syserr (43, "Null pointer indirection.\n"
                "It's an internal error. Ask developers.");
}

void syserr50 (void)
{
    syserr (50, "Integer overflow detected.\n"
                "It's an internal error. Ask developers.");
}

/* -----------------------------------------------------------------------
    Internal trace counters
----------------------------------------------------------------------- */

int pcount [PCOUNT];

/* -----------------------------------------------------------------------
    Name:       pcinit

    Purpose:    Initilaize values of counters.

    Usage:      pcinit()
----------------------------------------------------------------------- */

void pcinit (void)
{   
    int i;

    /* Initialize counters */

    for (i = 0; i < PCOUNT; i++)
    {
        pcount [i] = 0;
    }
}

/* -----------------------------------------------------------------------
    Name:       ptrace

    Purpose:    Increase trace counter.

    Usage:      ptrace (n)

    Returns:    Counter value.

    Arguments:
       n      - Counter
----------------------------------------------------------------------- */

#ifndef ptrace  /* Might be implemented as a macro in fexact.h */

int ptrace (int n)
{   
    /* Increase counter */

    return (++pcount [n]);
}

#endif

/* -----------------------------------------------------------------------
    Name:       pctype

    Purpose:    Type profiling counters to console.

    Usage:      pctype()
----------------------------------------------------------------------- */

void pctype (void)
{   
    int i;

    /* Type carriage return */

    fprintf (stderr, "\r");

    /* Type counters */

    for (i = 0; i < PCOUNT; i++)
    {
        if (pcount [i]) fprintf (stderr, "%d:%d ", i, pcount [i]);
    }
}
