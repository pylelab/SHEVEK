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

 

 
/** Standard Definitions **/

#define NO	0
#define YES	1


 /** ChiSquare and estimating probability routines -- from numerical recipies **/
 
#define TINY 1.0e-30    /*A small number.*/
#define ITMAX 1000		/* previously 100 */
#define EPS 3.0e-7
#define FPMIN 1.0e-30


/** Cochran Conditions Parameters -- determine how P value is calculated **/

#define PERCENT		0.80
#define EXPECTED	5.0
#define MINIMUM		1.0

 
/** Parameter limiting the number of different characters possible in a given position **/

#define POSITIONS	40


/** Statistical calculations are not always sensible, given certain data sets **/

#define MEANINGLESS -1     /*If statistical calculations are not possible, value set to -1*/
							

/** Used when calling chisquare or arcalc functions -- resets matrix orgin to 1,1 **
/** Eventually, this can be removed by redefining matrix which uses this		**/

#define ORIGIN 1

/** Multiplier for distributions **/

#define PMULT 4
#define VMULT 1

/** USED in threshold.c as definition of number of columns of data output in data files **/
 
#define COLS  8			/*NUMBER OF COL IN FILE*/

#define col1  5			/*COLUMN FOR CRAMER'S V SCORE*/
#define col2  7			/*COLUMN FOR - LOG P SCORE*/
#define col3  6			/*COLUMN FOR DF			*/
#define col4  2			/*previous COLUMN for X2, set to GROUP NUMBER*/
#define col5  3         /*previuous COLUMN for x2-approx P; set to Sensitivity*/

#define ALLOFFSET 3		/*number of textlines at end of alldata file -- see chi_analysis*/

/** Parameters used in misalignment functions **/


#define PCOLS  6				/*Number of columns in predictions.txt file -- see threshold.c */
#define MINIMUMSEQ 2
						 	/*Number of times sequence appears to be considered misaligned*/
#define ARTHRESH 1.0		/*Threshold for AR cellvalues -- distinguish misaligned from not*/
#define NOCHAR '#'			/*Character default if no suggested alternatives are possible */


/** Variable as flag -- when wish to indicate that a variable has meaningless or unused value **/

#define NEGONE  -1


/** Parameter for error and min on logP and V **/

#define VERROR 0.05			/*currently unused parameter*/
#define VMINIMUM 0.40
#define VMAXIMUM 1.0

#define PERROR 0.5			/*determined by area calcuations for bad predicts using 3D plot*/
#define PMINIMUM 2.0		/*see also BC calculated by Bonforroni correction: see distributionanalyzer.c*/
#define RESDEFAULT 0.1		/*resolution for progressive threshing*/
#define PSTABILITY 5		/* five units of RESDEFAULT*/

#define DFDIMERR 1
#define DFMAX	100.0

#define GAPAMNT 0.90

#define SDEFAULT 1.5



/** Flag for failure of FEXACT algorithm **/
/*(since probability is always 0 to 1)*/

#define EFAIL	2



/* DEFINTIONS for redundanteliminator.c */

#define MAXWIDTH 200
#define MAXTOKENS 10000				/*defines total number of redudant predictions possible*/



/* DF Imbalance Parameter */
#define PHYLO 3.0


/*************************/