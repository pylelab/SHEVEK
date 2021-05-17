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
/*	Certain algorithms found within this file may be derivatives*/
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
/************   CHIPANG function     ************************/
/**************************************************************/
  
float chipang(int **nn, int ni, int nj, float *chisq, float *df, double *prob, double *cramrv, float *ccc)

/*calculates chi, V, df, and sensitivity*/
/*derived from num recipie base code*/
/*Given a two-dimensional contingency table in the form of an integer array nn[1..ni][1..nj]*/

{

/*double gammq(double a, double x);  */
int nnj,nni,j,i,minij;
float sum=0.0,*sumi,*sumj;

double expctd, temp;

int ones = 0;

extern float CRV_CUTOFF;
extern int NUM_ROW;

float **xmat;
float temp2;
float xmatmax = 0;

float svalue;

double isum;

xmat = matrix(1, ni, 1, nj);     /*allocates memory for an int matrix */ 


sumi=vector(1,ni);
sumj=vector(1,nj);
nni=ni; 							/*Number of rows*/
nnj=nj; 							/*and columns.*/


/***********************************************************************/
/* Test Matrix for zero rows or columns - eliminate from df calculation -- redundant? check? gap_rectifier*/
/***********************************************************************/


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

/********************************************************************************************/
/*TEST matrix dimensions after zero row/column removal - must be min. 2 for rows and columns*/
/********************************************************************************************/
 

if ((nni < 2) || (nnj < 2))	{		/** Chisquare statistical analysis not possible **/
	*prob = MEANINGLESS;
	return 0;
}

/***************************************************************/
/*TEST OF GAPS in the positions -- N versus SUM for cramer's V */
/***************************************************************/


if ( sum <= ( ((float) NUM_ROW) * ((float) GAPAMNT) ))	{ 	/*i.e. this matrix represents less then %GAPANMT */
		*cramrv = 0.0;			  
		*prob = MEANINGLESS;	  
}								  
 

/****************************************/
/*TEST OF DF IMBALANCE in the positions */
/****************************************/


if ((nni >= ((float)PHYLO)*nnj) || (nnj >= ((float)PHYLO)*nni))	{
	*prob = MEANINGLESS;
	return 0;
	 
}



/*****************************************************************************************/
/*TEST matrix for singles - one for both column and row sum - part of sensitivity measure*/
/*****************************************************************************************/


for (i=1;i<=ni;i++) { 					 
	if (sumi[i] == 1.0) {
	for (j=1;j<=nj;j++) {
			if ( sumj[j] == 1.0 && nn[i][j] == 1) {  
				++ones;
				break;
			}
		}
	}
}


 
/**************************************/
/* parameter/statistical calculations */
/**************************************/


*df= (float) (nni*nnj-nni-nnj+1); 					/*Corrected number of degrees of freedom.*/
*chisq=0.0;

isum = 1/sum;		/*reduce number of division calculations*/

for (i=1;i<=ni;i++) { 					/*Do the chi-square sum.*/
	for (j=1;j<=nj;j++) {
		expctd=sumj[j]*sumi[i] * isum;
		temp=nn[i][j]-expctd;
temp2 = (float) (temp*temp/(expctd+TINY));
		xmat[i][j] = temp2;	
		
		/*printf("chi %f\n", xmat[i][j]); */	
		
								/***allocate matrix of chi values for each cell*/
		*chisq += temp2; 	
 
	}							/* eliminated row or column will*/
}								/*not contribute to the sum.*/



/* *prob=gammq(0.5*(*df),0.5*(*chisq)); */		/*Chi-square probability function.*/
/* DO NOT calc x2 approx prob for now */

*prob = 1;

minij = nni < nnj ? nni-1 : nnj-1;
*cramrv=sqrt(*chisq/(sum*minij));			
*ccc=(float) sqrt(*chisq/(*chisq+sum));



/************************************/
/*** Calculate Matrix Sensitivity ***/
/************************************/

	/*find Sij and Smax*/

if (*cramrv >= VMINIMUM) {   /*calc sensitivity for only higher V value scores */

	for (i=1;i<=ni;i++) { 				
		for (j=1;j<=nj;j++) {
			if (nn[i][j] > 0)	{

				xmat[i][j] = xmat[i][j] / ((*chisq)  * nn[i][j]);		/***set each cell equal to sensitivity score*/
			
				/*printf("score %f\n", xmat[i][j]);*/

				if (xmatmax < xmat[i][j]) xmatmax = xmat[i][j];
			
				/* printf("best %f\n",xmatmax); */
			}

		}
	}

}

svalue = (float) (((*cramrv) * xmatmax) + ones);  /*relative to V plus number of singles in matrix*/


/*********************/
/* Deallocate memory */
/*********************/


free_vector(sumj,1,nj);
free_vector(sumi,1,ni);

return svalue;

}


/** OLD CODE **/




/**************************************************/
/********FUNCTION: gammq               ************/
/*************************************************/

/*
double gammq(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) 	/*nrerror("Invalid arguments in routine gammq");*/
/*			return MEANINGLESS;
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

*/
/**************************************************/
/********FUNCTION: gcf                 ************/
/*************************************************/
/*

void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

*/
/**************************************************/
/********FUNCTION: gser                ************/
/**************************************************/
/*
void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
*/
/**************************************************/
/********FUNCTION: gammln              ************/
/**************************************************/
/*
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


*/


/**************************************************************/
/************   CHISQUARE function     ************************/
/**************************************************************/
  
/*  int chisquare(int **nn, int ni, int nj, float *chisq, float *df, double *prob, double *cramrv, float *ccc)

/*Given a two-dimensional contingency table in the form of an integer array nn[1..ni][1..nj],
this routine returns the chi-square chisq, the number of degrees of freedom df, the signi_cance
level prob (small values indicating a signi_cant association), and two measures of association,
Cramer's V (cramrv) and the contingency coe_cient C (ccc).  */

/*Since you pass the address, and fill the value, then of course, the address you passed, will have a value in it!*/
/*this is how multiple values are returned! */

/*

{
double gammq(double a, double x);
int nnj,nni,j,i,minij;
float sum=0.0,*sumi,*sumj;

double expctd, temp;

extern float CRV_CUTOFF;
extern int NUM_ROW;

sumi=vector(1,ni);
sumj=vector(1,nj);
nni=ni; 							/*Number of rows*/
/*  nnj=nj; 							/*and columns.*/
/*
for (i=1;i<=ni;i++) { 					/*Get the row totals.*/
/*	sumi[i]=0.0;
	for (j=1;j<=nj;j++) {
		sumi[i] += nn[i][j];
		sum += nn[i][j];
	}
	if (sumi[i] == 0.0) --nni; 			/*Eliminate any zero rows by reducing the num*/
/*}


for (j=1;j<=nj;j++) { 					/*Get the column totals.*/
/*	sumj[j]=0.0;
	for (i=1;i<=ni;i++) sumj[j] += nn[i][j];
	if (sumj[j] == 0.0) --nnj; 			/*Eliminate any zero columns.*/
/*}




/**********************************************************/
/*** 2nd Test of Matrix: must be 2x2 after zero removal ***/
/**********************************************************/
/*
if ((nni < 2) || (nnj < 2))	{		/** Chisquare statistical analysis not possible **/
/*	*prob = MEANINGLESS;
	return 0;
}


/**********************************************************/
/*****************  End of 2nd Test ***********************/
/**********************************************************/


/*
*df= (float) (nni*nnj-nni-nnj+1); 					/*Corrected number of degrees of freedom.*/
/*  *chisq=0.0;

for (i=1;i<=ni;i++) { 					/*Do the chi-square sum.*/
/*	for (j=1;j<=nj;j++) {
		expctd=sumj[j]*sumi[i]/sum;
		temp=nn[i][j]-expctd;
		*chisq += (float) (temp*temp/(expctd+TINY)); 	/*Here TINY guarantees that any*/ /*elim. Div. By zero */
/*	}							/* eliminated row or column will*/
/*}								/*not contribute to the sum.*/

/*
*prob=gammq(0.5*(*df),0.5*(*chisq)); 		/*Chi-square probability function.*/
/* minij = nni < nnj ? nni-1 : nnj-1;
*cramrv=sqrt(*chisq/(sum*minij));			
*ccc=(float) sqrt(*chisq/(*chisq+sum));



/********* TEST OF N versus SUM for cramer's V **********/
/*
if ( sum < (NUM_ROW/2.0))	{ 	/*i.e. this matrix represents less then 50% of sequences */
/*		*cramrv = 0.0;			  
		*prob = MEANINGLESS;	  
}								  

/*********************************************************/

/*
free_vector(sumj,1,nj);
free_vector(sumi,1,ni);

return 0;
}

*/
