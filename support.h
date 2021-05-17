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
    Heder file for support functions.
----------------------------------------------------------------------- */

extern void* sysalloc(int nitems, int size);
extern void  sysfree (void *p);
extern void  syschk  (void);

extern void syserr   (int icode, const char *mes);
extern void syserr40 (void);
extern void syserr41 (void);
extern void syserr42 (void);
extern void syserr43 (void);
extern void syserr50 (void);

#define PCOUNT 100              /* Number of profiling counters */
extern int pcount [PCOUNT];         /* Profiling counters */

extern void pcinit (void);
extern int  ptrace (int n);
extern void pctype (void);



#define ptrace(n) (++pcount[n])

/* -----------------------------------------------------------------------
    Name:       RAND

    Purpose:    Generate random value in the range [0..1)

    Usage:      RAND

    Result:     The random value.
----------------------------------------------------------------------- */

/* #define RAND ((float) rand() / (float) RAND_MAX)  */

/* #define RAND PANGRAND(); */

