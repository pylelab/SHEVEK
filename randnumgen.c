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


/*ORIGINALLY TAKEN FROM:

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
/*																	   */
/*		Portable Random Number Generators							   */
/*																	   */
/*		Gerald P. Dwyer Jr. and K.B. Williams						   */
/*		Federal Reserve Bank of Atlanta								   */
/*		Working Paper 99-14										       */
/*		October 1999												   */
/*																       */
/*			NOTE: PORTIONS OF THIS CODE IN THIS FILE have been taken   */
/*			from the above paper									   */
/* 			the above paper.										   */
/* 																	   */
/* 																	   */
/*																	   */
/***********************************************************************/

/* CONTACT THE ABOVE AUTHORS AT:

Gerald P. Dwyer Jr., 
Research Department, 
Federal Reserve Bank of Atlanta, 
104 Marietta Street, NW, 
Atlanta, Georgia 30303-2713, 
404/614-7095, 
404/521-8810 (fax), 
gerald.p.dwyer@atl.frb.org 
(or dwyerg@clemson.edu); 

K.B. Williams, 
412 Magnolia Avenue, 
Melbourne Beach, Florida 32951-2000
kbwms@aol.com. 

*/


/* THIS PROGRAM WAS THEN OPTIMIZED AS A WORK FOR HIRE PROJECT by RAUL Shakirov */
 



/* THE FOLLOWING KNOWN BUGS WERE ELIMINATED **********************************************/
/*
  1. pangrold.c line 83

   static const long r2 = 2147483587 - 44095 *  (2147483647 / 65670);   mod2 - mult2 * q1 

   Problem:  Mistyping. Must be:

   static const long r2 = 2147483587 - 44095 *  (2147483647 / 44095);   mod2 - mult2 * q2 

   Status:   Fixed in both original and optimized code.


  2. pangrold.c line 83

   randtemp1 =   ((long double) RandComb())  /  ( (long double) (mod1-1) );

   Problem:  RandComb Provides for value in the rande [0..,mod1-2].
             If one need long double value in the range [0..1]
             one must divide by mod2-2:

   randtemp1 =   ((long double) RandComb())  /  ( (long double) (mod1-2) );

   Status:   Fixed in both original and optimized code.
*/
/******************************************************************************************/

/*** NOTE: **********************************************************************************/
/*
	1) The Estimated Cycle Length is 2.3 x 10^18 
	
	2) The Start has been set to 1 for each generator; this "seed" could be changed each time
*/
/*********************************************************************************************/

/* -----------------------------------------------------------------------


    Combination random generator with approximate factoring.

    Language:     ANSI C
    Arithmetic:   32-bit
----------------------------------------------------------------------- */

#include <math.h>

/*
   combination multiplicative random number generator
   takes difference between two random numbers:

   rand1 = (rand1 * mult1) % mod1   (range 1..mod1-1)
   rand2 = (rand2 * mult1) % mod2   (range 1..mod2-1)
   rand  = (rand1 - rand2) % mod1   (range 0..mod1-2)

   rand1 and rand2 can be initialized outside of file
*/

#define mod1    (2147483647L)
#define mult1   (65670L)
#define q1      (mod1 / mult1)
#define r1      (mod1 - mult1 * q1)

#define mod2    (2147483587L)
#define mult2   (44095L)
#define q2      (mod2 / mult2)
#define r2      (mod2 - mult2 * q2)

#define start1  (1)
#define start2  (1)

#define LONG_BITS   (8 * sizeof (long) - 1)



/* -----------------------------------------------------------------------
    Optimized generator for 32-bit arithmetic, no conditional branches.
----------------------------------------------------------------------- */

long double pangrand2 (void)
{
    static long rand1 = (start1 * mult1) % mod1;
    static long rand2 = (start2 * mult2) % mod2;

    long temp1 = rand1;     long temp2 = rand2;
    long temp1s = q1;       long temp2s = q2;

    /*
        Compute return value
    */

    long temp = temp1 - temp2;
    /*
        if (temp < 0) { temp += mod1; }
    */
    temp += ((temp >> LONG_BITS) & (mod1 - 1));

    /*
        Compute new random values
    */

    temp1 /= q1;            temp2 /= q2;

    temp1s *= temp1;        temp2s *= temp2;
    temp1s -= rand1;        temp2s -= rand2;
    temp1s *= -mult1;       temp2s *= -mult2;
    temp1 *= r1;            temp2 *= r2;
    temp1s -= temp1;        temp2s -= temp2;
    /*
        if (temp1s < 0) { temp1s += mod1; }
        if (temp2s < 0) { temp2s += mod2; }
    */
    temp1 = temp1s;         temp2 = temp2s;
    temp1s >>= LONG_BITS;   temp2s >>= LONG_BITS;
    temp1s &= mod1;         temp2s &= mod2;
    temp1s += temp1;        temp2s += temp2;

    /*
        Store new random values
    */

    rand1 = temp1s;         rand2 = temp2s;

    /*
        Convert to long double in the range [0..1]
    */

    return ((long double) temp) / (mod1 - 2);
}


 















/** ORIGINAL OLD PANGRAND FUNCTION ***/



/***********************************************************************/
/***********************************************************************/


/* RandComb.cpp
combination multiplicative random number generator
takes difference between two random numbers
rand1 and rand2 can be initialized outside of file */
/*

static const long mod1 = 2147483647 ;
static const long mult1 = 65670 ;
static const long q1 = 2147483647 / 65670;  /* mod1/ mult1 ;*/
/*static const long r1 = 2147483647-65670 * (2147483647 / 65670);  /*mod1 - mult1 * q1 ;*/

/*
static const long mod2 = 2147483587 ;
static const long mult2 = 44095 ;
static const long q2 = 2147483587/44095; /* mod2 / mult2 ; */
/*static const long r2 = 2147483587 - 44095 *  (2147483647 / 65670); /*mod2 - mult2 * q1 ;*/
/*

long rand1=1, rand2=1 ;



/**************************************************************/
/************   PANGRAND function					 **********/
/**************************************************************/
/*long double pangrand(void)
{
	long double randtemp1;
	
	/*long double randtemp2;*/

/*	long RandComb(void);

	randtemp1 =   ((long double) RandComb())  /  ( (long double) (mod1-1) );

	/**randtemp2 = randtemp1 /  ( (long double) (mod2-1) );  **/

/*	return randtemp1;  /*  2*randtemp2;  */
/*
}

/**************************************************************/
/************   GENRRAND function					 **********/
/**************************************************************/


/* generate the next value in sequence from generator using approximate factoring */
/*
long GenrRand(long rand, long modulus, long mult, long q, long r)
{

	long temp ;
	
	temp = rand / q ;
	temp = mult * (rand - temp * q) - temp * r ;
	
	if (temp < 0)	temp += modulus ;
	
	return temp ;
}




/**************************************************************/
/************   RANDCOMB function					 **********/
/**************************************************************/

/* get a random number */
/*
long RandComb(void)
{
	long RandNum ;
	long GenrRand(long rand, long modulus, long mult, long q, long r);


	rand1 = GenrRand(rand1, mod1, mult1, q1, r1) ;
	rand2 = GenrRand(rand2, mod2, mult2, q2, r2) ;
	RandNum = rand1 - rand2 ;
	
	if (RandNum < 0)	RandNum += mod1-1 ;
	
	return RandNum ;
}*/