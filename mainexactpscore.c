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

 

/******************************************/
/******************************************/
/* Estimated Exact Probability Function  **/
/* via Monte Carlo Simulation            **/
/*										 **/
/* For use when exact (FEXACT function)  **/
/* Fails; and Cochran Conditions are not **/
/* met									 **/
/******************************************/
/******************************************/


#include <math.h>
#include <stdio.h>

#include "nrutilp.h" 
#include "rcont2p.h" 
#include "shevek.h"
#include "definitions.h"



/**************************************************************/
/**************************************************************/
/************   ESTEXACT function       ***********************/
/**************************************************************/
/**************************************************************/
/**************************************************************/
/************   ESTEXACT function       ***********************/
/**************************************************************/
/**************************************************************/
/**************************************************************/
/************   ESTEXACT function       ***********************/
/**************************************************************/
/**************************************************************/


/* returns p, no matter how calculated */

double estExact(int *ROWmatrix, double *dmatrix, double PRECISION, int numrows, int numcols, int *Echeck)
{

int rcont2s(int nrow, int ncol, int *_nrowt,int *_ncolt,int *_matrix);
double rcont2f (void);
/* int	fexactmain (int nrow, int ncol, double *table, double *Eprob);*/
 
int minprecision = 500;    /* tested 10,000; 2000 works as well; 1/2000 ~ 5e-04 */

/*double *Eprob;  /* exact probability as calculated by FEXACT algorithm method */
/* double Eprobmem;	
int Eflag = 0;		/* flag determining of Eprob value is good = 0 */
/*double Etemp;*/

int *RANDmatrix;
int *rowtot; 
int *coltot; 
float *expctd;

double probability;

double temp;
double chisq;
double chisqRAND;
double tprob;				/*VARIABLE IS UNUSED FOR CALCULATIONS*/
int tally;
float ftally;
int numrows2,numcols2,i,j,k,m,n,r, error;
float sum=0.0;

/***********************/
/*allocate memory      */
/***********************/
 
/* Eprob = &Eprobmem; */
RANDmatrix=ivector(0,(numrows*numcols -1));
rowtot=ivector(0,(numrows-1));
coltot=ivector(0,(numcols-1));
 
expctd=vector(0,(numrows*numcols -1));


/***********************/
/*clears memory = 0    */
/***********************/


for (r=0;r < numrows*numcols; r++) {
	*(RANDmatrix + r) = 0;
	*(expctd + r) = 0;
}


/***********************/
/*calc row and col info*/
/***********************/
 
numrows2=numrows; 							/*Number of rows*/
numcols2=numcols; 							/*and columns.*/

for (i=0;i< numrows;i++) { 					/*Get the row totals.*/
	rowtot[i]=0;							/*clears memory*/
	for (j=0;j< numcols;j++) {
		rowtot[i] += *(ROWmatrix +(i*numcols + j));
		sum += *(ROWmatrix +(i*numcols + j));
	}
	if (rowtot[i] == 0) --numrows2; 			/*Eliminate any zero rows by reducing the num*/
}


for (j=0;j< numcols;j++) { 					/*Get the column totals.*/
	coltot[j]=0;							/*clears memory*/
	for (i=0;i< numrows;i++) coltot[j] += *(ROWmatrix +(i*numcols + j)); 
	if (coltot[j] == 0) --numcols2; 			/*Eliminate any zero columns.*/
}


/**********************************************************/
/*** Test of Matrix: must be 2x2 after zero removal *******/
/**********************************************************/

if ((numcols2 < 2) || (numrows2 < 2))	{		
	/** Chisquare statistical analysis not possible **/
	return MEANINGLESS;
}


/************************************/
/* Calculate chi of original Rmatrix*/
/************************************/

chisq=0.0;

for (i=0;i< numrows;i++) { 					/*Do the chi-square sum.*/
	for (j=0;j< numcols;j++) {
		*(expctd +(i*numcols + j)) = (coltot[j]*rowtot[i])/sum;
		temp= *(ROWmatrix +(i*numcols + j)) - *(expctd +(i*numcols + j));
		chisq += temp*temp  /  (*(expctd +(i*numcols + j)) + TINY); 	
								/*Here TINY guarantees that any*/ /*elim. Div. By zero */
	}							/* eliminated row or column will*/
}								/*not contribute to the sum.*/




/****************************************************/
/* MAIN CYCLE										*/
/*									   				*/
/* Three parts:										
  
   Rationale: 

  Montecarlo simulation requires immense computational power;
  however, since the precision of the calculated p value
  is related to 1/simulations, then for high p values (close to 1)
  it stands to reason that not many simulations are necessary,
  since ultimately the log of the p value is taken and it is
  the order of magnitude that is most important.

  Thus, the montecarlo simulation is broken into two; if the
  p value calculated is high, then the simulation stops at 2000
  simulant tables.

  If however, the p value is low, then it must be 
  further simulated.  

  The reason a small MC simulation is first attempted is that
  such a small simulation will be much faster than an actual exact
  calculation -- in other words, the exact method is far superior
  to simulation in terms of speed, as the p value becomes lower.

  THUS: 1) mini precision MC simulation
	 
		2) full MC simulation
/*													*/
/*													*/ 
/****************************************************/


/** added to jump over scoring -- see shevek.c file and openfile function **/
/** programer's hack -- MUST see openfile() function for other required skip**/

if (PRECISION == 42)	{	 
	minprecision = 16;			/* jump over MC simulation */
}	

if (PRECISION == 43)	{	 
	minprecision = 16;			/* jump over FEXACT calculation */
}	

/**************************************************************************/

tally = 0;
error = 0;
error = rcont2s(numrows,numcols,rowtot,coltot,RANDmatrix);  /* test: is matrix good? */

if (error == 0) {

	/************/
	/** PART I **/
	/************/



	for (k=1; k<=minprecision; k++)	{


		/*********************************************************/
		/*** generate random matrix using marginal frequencies ***/
		/*********************************************************/
 
		tprob = rcont2f();

		/**************************************************/
		/*** determine chisquare for that random matrix ***/
		/**************************************************/

		chisqRAND = 0.0;

		for (m=0;m< numrows;m++) { 													/*Do the chi-square sum.*/
			for (n=0;n< numcols;n++) {
				temp= *(RANDmatrix +(m*numcols + n)) - *(expctd +(m*numcols + n));
				chisqRAND += temp*temp  /  (*(expctd +(m*numcols + n)) + TINY); 		/*Here TINY guarantees that any*/ 
																					/*elim. Div. By zero */
			}																		/* eliminated row or column will*/
		}																			/*not contribute to the sum.*/


		/**************************************************/
		/* COMPARE original matrix to random matrix       */
		/**************************************************/


		if (chisqRAND >= chisq)
			tally += 1;
		
	}





		/*************
	
***REMOVED UNTIL FURTHER NOTICE **/
	
		/** PART II **/
		/*************/ 
	*Echeck;
	*dmatrix;
	/*printf("tally %d", tally);*/


/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/
/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/
/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/
/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/
/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/
/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/ 
/* ? QUESTION: should Echeck be a criteria to enter this loop??? **/
	/* determined by analysis of error 63 */

/*
	if ((tally < 10) && (minprecision != 16))	{	
		/* i.e. p is not > than 10/miniprecision i.e. "p is low" */
		/* and programer does not wish to skip FEXACT calculation */
/*
		Eflag = fexactmain(numrows, numcols, dmatrix, Eprob);
	
		if (Eflag == 0)		{		/*calcuation is sucessful*/
/*			Etemp = *Eprob * (-1);
			return Etemp;			/*send back FEXACT probability; -1 acts as counter */
/*		}

		else	{					/*FAILURE FAILURE FAILURE */
						 			/*signficant since all very low P need to be calculated*/
									/*with same possible maximum*/
/*
					printf("\n    WARNING WARNING WARNING WARNING\n");
					printf("	FAILURE OF EXACT P CALCULATION:\n");
					printf("	DATA SET MOSTLIKELY LARGE.     \n");
					printf("	User should Limit logP scores. \n");
					printf("    (Max should be 1/PRECISION.)   \n\n");
					printf("	NOTE: Calculations must now be \n");
					printf("		  determined by simulation.\n");
					printf("		  MAY TAKE DAYS OR MONTHS. \n");
					printf("		  To END, hit CNTR-C.	   \n");

					*Echeck = 1;
				
		}

	}

 


	/**************/
	/** PART III **/
	/**************/ 

		 
	while ((tally<10) && (k<=PRECISION))	{

		/*********************************************************/
		/*** generate random matrix using marginal frequencies ***/
		/*********************************************************/
 
		tprob = rcont2f(); 

		/**************************************************/
		/*** determine chisquare for that random matrix ***/
		/**************************************************/

		chisqRAND = 0.0;

		for (m=0;m< numrows;m++) { 													/*Do the chi-square sum.*/
			for (n=0;n< numcols;n++) {
				temp= *(RANDmatrix +(m*numcols + n)) - *(expctd +(m*numcols + n));
				chisqRAND += temp*temp  /  (*(expctd +(m*numcols + n)) + TINY); 		/*Here TINY guarantees that any*/ 
																					/*elim. Div. By zero */
			}																		/* eliminated row or column will*/
		}																			/*not contribute to the sum.*/


		/**************************************************/
		/* COMPARE original matrix to random matrix       */
		/**************************************************/


		if (chisqRAND >= chisq)
			tally += 1;

		k = k+1;

	}


}


else return MEANINGLESS;


 
/*******************/
/*END OF MAIN CYCLE*/
/*******************/


if (tally != 0) {
	ftally = (float) tally;
	probability = ftally/ (k-1);
}
else
	probability = 1.0/PRECISION;


free_ivector(coltot,0,(numcols-1));
free_ivector(rowtot,0,(numrows-1));
free_ivector(RANDmatrix,0,(numrows*numcols -1));
free_vector(expctd,0, (numrows*numcols -1));
 

return probability;


}





 