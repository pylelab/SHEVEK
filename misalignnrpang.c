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



/****************************************************************/
/****************************************************************/
/****************************************************************/
/*																*/
/*	The algorithms found within this file may be derivatives	*/
/*  of source code obtained from the book:						*/
/*  "Numerical Recipes in C: The Art of Scientific Computing"	*/
/*  published by Cambridge University Press.					*/
/*																*/
/****************************************************************/
/****************************************************************/
/****************************************************************/
 
#include <math.h>

#include "nrutilp.h"
#include "definitions.h" 


/**************************************************************/
/************   ARCALC function     ***************************/
/**************************************************************/
 

void arcalc(int **nn, int ni, int nj, float *chisq, float *df, double *prob, float *cramrv, float *ccc, float **ar)

/*Given a two-dimensional contingency table in the form of an integer array nn[1..ni][1..nj],
this routine returns the chi-square chisq, the number of degrees of freedom df, the signi_cance
level prob (small values indicating a signi_cant association), and two measures of association,
Cramer's V (cramrv) and the contingency coe_cient C (ccc).  */

/*Since you pass the address, and fill the value, then of course, the address you passed, will have a value in it!*/
/*this is how multiple values are returned! */

{
/*float gammq(float a, float x);*/
int nnj,nni,j,i;
float sum=0.0,expctd,*sumi,*sumj,temp;

double temp2;			/****for AR calc***/

ccc;
cramrv;
prob;
df;

sumi=vector(1,ni);
sumj=vector(1,nj);
nni=ni; 							/*Number of rows*/
nnj=nj; 							/*and columns.*/


for (i=1;i<=ni;i++) { 					/*Get the row totals.*/
	sumi[i]=0.0;
		for (j=1;j<=nj;j++) {
			sumi[i] += nn[i][j];
			sum += nn[i][j];
		}
	if (sumi[i] == 0.0) --nni; 			/*Eliminate any zero rows by reducing the num*/
}


for (j=1;j<=nj;j++) { 					/*Get the column totals.*/
	sumj[j]=0.0;
	for (i=1;i<=ni;i++) sumj[j] += nn[i][j];
	if (sumj[j] == 0.0) --nnj; 			/*Eliminate any zero columns.*/
}

 /* *df=nni*nnj-nni-nnj+1;  */ 			/*Corrected number of degrees of freedom.*/
*chisq=0.0;

for (i=1;i<=ni;i++) { 					/*Do the chi-square sum.*/
	for (j=1;j<=nj;j++) {
	expctd=sumj[j]*sumi[i]/sum;
	temp=nn[i][j]-expctd;
	/*  *chisq += temp*temp/(expctd+TINY); 	*/  /*Here TINY guarantees that any*/ /*elim. Div. By zero */

	temp2 = expctd * (1 - sumi[i]/sum) * (1 - sumj[j]/sum);

	ar[i][j] = (float) ( (temp) / sqrt(temp2)  );     /***AR calc***/

		/*** printf("\n%f %f\n", temp2, ar[i][j]);	***/

	}							/* eliminated row or column will*/
}	

							/*not contribute to the sum.*/

	/* *prob=gammq(0.5*(*df),0.5*(*chisq)); 	*/	/*Chi-square probability function.*/
	/* minij = nni < nnj ? nni-1 : nnj-1;       */
	/* *cramrv=sqrt(*chisq/(sum*minij));		*/
	/* *ccc=sqrt(*chisq/(*chisq+sum));			*/

free_vector(sumj,1,nj);

free_vector(sumi,1,ni);
}
