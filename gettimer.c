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


#include <time.h>       /* Header file for standard time functions */

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

/* -----------------------------------------------------------------------
    Name:       GetFloatTimer

    Purpose:    Get timer as a second interval between current time
                and given time value.

    Usage:      float timer = GetFloatTimer (0);
                ...
                Some time-consuming process
                ...
                timer = GetFloatTimer (timer);

    Arguments:
      timer   - Time value to be substracted from current time.
----------------------------------------------------------------------- */

float GetFloatTimer (float timer)
{
    /*
        To get correct time difference we need double computation,
        because clock() value may be too great for float presision!
     */

    float time = (float) (((double) clock() / CLK_TCK) - timer);

    /*  Correct time after midnight  */

    if (time < -1.) time += (3600 * 24);

    /*  Correct impossible precision error. */

    if (time <  0.) time = 0.;

    return time;
}
