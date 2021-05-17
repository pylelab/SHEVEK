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




#include <stdio.h>      /* Header file for standard io functions */
#include <stdlib.h>     /* Header file for standard functions */
#include "rcont2p.h"     /* Header file for rcont2 functions */

/* -----------------------------------------------------------------------
    Name:       prcerr

    Purpose:    Print an error message to console and exit.

    Usage:      prcerr (icode, mes)

    Arguments:
      icode   - Int code for the error message.                 (Input)
      mes     - Character string containing the error message.  (Input)
----------------------------------------------------------------------- */

void prcerr (int icode, const char *mes)
{
    fprintf (stderr, "\nRCONT2 ERROR: %d %s\n", icode, mes);
                                /* Print diagnostic to console */
    exit (2);                   /* Close all files and exit with code 1 */
}

/* -----------------------------------------------------------------------
    Name:       prcerrN

    Purpose:    Print error message N to console and exit.

    Usage:      prcerrN()
----------------------------------------------------------------------- */

void prcerr1 (void)
{
    prcerr (1, "Number of rows must be greater than 1.");
}

void prcerr2 (void)
{
    prcerr (2, "Number of columns must be greater than 1.");
}

void prcerr3 (void)
{
    prcerr (3, "Row totals must be greater than 0.");
}

void prcerr4 (void)
{
    prcerr (4, "Column totals must be greater than 0.");
}

void prcerr5 (void)
{
    prcerr (5, "Sample size exceeds 5000.");
}

void prcerr6 (void)
{
    prcerr (6, "Sum of row and column totals must be the same.");
}
