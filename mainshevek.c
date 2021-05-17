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
 



 
/** NOTES, BUGS, CONSIDERATIONS **/



/* NOTE: redundant eliminator compares strings, not numbers, 
therefore is not in 'numerical order' -- but who cares? */

/* NEED TO PUT IN PROGRAMMER"S LOGIC BUGS -- just a check everywhere */

/* 2-12-02 fixed bug in eliminator2 -- groupnumber < 0! 
made it so that need to check for neg1 first */


/* STACK SIZE : set to 150MB committed; with 50MB reserve 
should let user know */


/** PUT INTO MANUAL:

B) if V is negative -2+ then more than one member of chain

C) if S is followed by 666 then "unstable prediction"**/


 
/*** Get rid of division where possible ***/

/** incorporate FEXACT **/
/** incorporate optimized estExact **/
/** build version which pulls in list -- so that can do first and second scan **/

/******************/

/* NOTES: 5-21-02

  Removed FEXACT: major problem: too buggy, with unknown logic path; cannot overcome
  coded limitations on dimensions/numbers -- i.e. the limitations appear to be most
  in the size/dimensions, not number of sequences. Should switch to MC Rescue methods.

  Added DF threshold and distribution calculator

  Will add today the DF IMBALANCE analyzer */


/***********************************************************************************/


/* NOTES 1-22-02
  Attempting to incorporate FEXACT algorithm; requires C and CPP files to work
  together; fails on linking: NOTE, need "/TC" command line to compile all as C files.
  (rather than C++ files)
*/


/***********************************************************************************/
/* FIXED FIXED FIXED FIXED FIXED FIXED FIXED FIXED FIXED FIXED FIXED FIXED *
/**** IMMEDIATE: GAP CHECK seems flawed -- results in any gap being bad??? ******/

/****** -- gap check -- problem issolated: In EstExact: b/c gap calc works in chipang, but not
in estExact -- not designed to handle zero rows and columns ***************/
/***************************************************************************

/* "size of distribution" -- i.e. number of elements check -- not by length, but by number of alldata points*/


/*  if more than one remains in chain elim, after all, eliminate all **/



/* STACK SIZE and recursive pathfinder function: because it is recursive,
memory allcation may be problematic; this is OS dependent; thus:
	FIRST: reduce stack requirements, by eliminating variables xup,yup,xdown,ydown
	and replace with loc variable; eliminate filesize variable with global one
	SECOND: increase stack size: in projects->settings->link->output	*/



/** version notes: ***/

/* 11-19-01 */

/* introduced optimized random number generator, by R.S. */
/* corrected rcont2 so that greater precision float --> double */



/*************************/
/*************************/

/**   PROGRAM BEGINS   ***/

/*************************/
/*************************/




/*** ORIENTATION NOTES *******************/
/*
	ON crosstab tables:


  matrix[p][q]
  matrix[i][j]
  matrix[row][column]

  position1 = count1 = n;p;ni;i = rows
  position2 = count2 = m,q,nj,j = columns

  NOTE: count1 and count2 start at 0; thus count1 = 1 means 2 rows;

*/
/******************************************/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "nrutilp.h"
#include "shevek.h"
#include "definitions.h"
 

double PRECISION; 
int prec;

float	NUM_COL;
int		NUM_ROW;
int		OFFSET;
 
int		mainseq;

int		START;
int		STOP;



/************************************************/
/************************************************/
/************   MAIN FUNCTION *******************/
/************************************************/
/************************************************/
/************************************************/
/************************************************/
/************   MAIN FUNCTION *******************/
/************************************************/
/************************************************/
/************************************************/
/************************************************/
/************   MAIN FUNCTION *******************/
/************************************************/
/************************************************/
/************************************************/
/************************************************/
/************   MAIN FUNCTION *******************/
/************************************************/
/************************************************/


int main ()
{
	void	intro();
	void	file_manager(void);
	char	*openfile();
	void	score_manager(char *input);
	void	positionrelater(char *input);
	void	screener(float *VThresh, float *PThresh, float *DFThresh);
	void	exitprogram();
	void	misalign_identifier(char *input, int offset, int numrows);


	char	*input;

	float VT,PT, DT;
	float *VThresh, *PThresh, *DFThresh;

	extern int prec;


	VThresh = &VT;
	PThresh = &PT;
	DFThresh = &DT;
	 

	intro();							/*INTRODUCTORY TEXT */

	file_manager();						/*handles/deletes any old data files that will conflict*/

	input = openfile();					/*OPENS THEN READS sequence alignment file*/
										/*returns pointer to sequence data*/

	positionrelater(input);				/*outputs correlation between position number and residue/nucleotide*/
										/*number of primary sequence under analysis*/


	score_manager(input);				/*SCORES sequence data positions*/
										/*outputs data files*/ 
										/*STANDARDIZES scores as necessary*/

	screener(VThresh,PThresh, DFThresh);/*SCREENING of Scores */
										/*1) Analyzes Distribution of Standardized scores*/
										/*2) Suggests lower min. thresholds*/
										/*3) Uses thresholds for preliminary list of predictions*/
										/*4) Eliminates intersecting interactions */
										/*5) Outputs predictions.txt file of predictions */ 
	
	
	misalign_identifier(input,OFFSET,NUM_ROW);	/*takes prediction list,and generates*/
												/*regenerates actual freq tables -- avoid having to keep them*/
												/*generates AR tables*/
												/*outputs misalign.txt file*/
								

	exitprogram();					


	free_cvector(input, 0, 0);		/*frees memory allocated for 	*/
										/*aligned sequences		*/

    return 0;
}


/*********************************************/
/************  INTRO  ************************/
/*********************************************/

void	intro()
{

   int  numread;
   int  numwrite;
   char buf[10];
   FILE *inputfile = NULL;

   /* opens input file */

	inputfile = fopen("titlepage", "r");		

	if (inputfile == NULL)	printf("\nTITLE PAGE NOT FOUND\n");
			   
	else {
	/*read FILE and writes to screen until end of file is reached*/

		while(  !feof(inputfile)   )
		{
		  numread =  fread(buf,sizeof(char),1,inputfile); 
		  numwrite = fwrite(buf,sizeof(char),1,stdout);
		} 
		
		fflush(inputfile);
	    fclose(inputfile);
	}
}
  
 

/*********************************************/
/************ openfile  **********************/
/*********************************************/

char	*openfile()
{

	char	*read_input(FILE *ifile);
	char	*input;
	char	filename[50];
 
	
	extern int	mainseq;
	extern double PRECISION;
	extern int prec;

	extern int START;
	extern int STOP;
 
	FILE * inputfile = NULL;

	int	found_file = 0;

	while (found_file == NO) {						/*Queries for input file*/
		printf("\nENTER alignment file name (include extension):");
		scanf("%s", filename);
		inputfile = fopen(filename, "r");
		if (inputfile == NULL)	 {
			printf("\nFAILURE TO OPEN OR LOCATE FILE!!!\n");
			found_file = NO;
		}
		else
		    found_file = YES;
	}


	printf("\n\nPlease enter the sequence number of the primary sequence:");
	scanf("%d",&mainseq);								/*query for mainseq parameter*/
 
	printf("\n\nEnter the INTEGER PRECISION of the P-values you will calculate(ex:4):");
	scanf("%d",&prec);									/*query for PRECISION parameter*/
	/* value of 42 -> "fast forwards" through scoring step */


	/* see score_manager function for use of these values below */
	printf("\n\nEnter the COLUMN POSITION start value (0-Total Positions):");
	scanf("%d", &START);

	printf("\n\nEnter the COLUMN POSITION STOP value (0-Total Positions; -1 = END):");
	scanf("%d", &STOP);

	
/* added to "jump" over score caculation step */

	if (prec == 42)	
		PRECISION = 42;

	if (prec == 43)
		PRECISION = 43;

	else
		PRECISION = (double) pow(10,prec); 

/**********************************************/




	input = read_input(inputfile);

	fclose(inputfile);

	return input;

}



/**************************************************/
/********FUNCTION: read_input *********************/
/** returns pointer to input **********************/
/**************************************************/


char	*read_input(FILE *ifile)
{
	char *cvector(long nl, long nh);

	void exitprogram();

	extern int		NUM_ROW;
	extern float	NUM_COL;
	extern int		OFFSET;

	char	*input;

	int	x;			/*int read*/
	char	c;			/*char read*/
	float	fmax_col;	/*max num of residues, float*/
	int	imax_col;		/*max num of residues, int*/
	int	seq = 0;			/*number of sequences = NUM_ROW*/

	int	count = 0;		/*total number of aa + label characters*/
	int	i = 0;		/*input index*/

	int	reading = 1;	/*flag for EOF*/

	long temp;

	NUM_COL = 0;
	NUM_ROW = 0;
	OFFSET = 0;



	printf("\nAllocating Memory for Input File.\n");

	/* LIMITED TO FASTA FORMAT INPUT DATA */
	/* LIMITED TO FASTA FORMAT INPUT DATA */
	/* LIMITED TO FASTA FORMAT INPUT DATA */
	
	x = fgetc(ifile);
	while (x != EOF) {
		if (x == '>')	{

			NUM_ROW += 1;		/* count the number of sequences */
			while (x != '\n')			/*exclude sequence names*/
				x = fgetc(ifile);
		
		}

		count += 1;
		x = fgetc(ifile);
	}


	fmax_col = (float) (count / NUM_ROW);
	imax_col = (int) fmax_col + 1;

	temp = NUM_ROW*imax_col;
	input = cvector(0,temp);
	/*allocate input matrix */
	/*NOTE: matrix will be larger than needed*/


	if (input == NULL) {
		printf("\nDynamic Memory Allocation FAILURE!!\n");
		exitprogram();
		exit(1);
	}
	else
		printf("Allocation Completed.\n");


	fseek(ifile, 0, SEEK_SET);	/*start reading from begining of file*/

	printf("Attempting to Read File. Please Wait.\n");
	x = fgetc(ifile);
	while (reading == YES)  {
		
		while (x != '>' && x != EOF )  /*scan past junk at begining of file*/
			x = fgetc(ifile);
		
		while (x != '\n'&& x != EOF )	/*scan pass sequence names*/
			x = fgetc(ifile);
		
		while ( x != '>' && x != EOF )   {			 /*See ASCII chart*/
			c = (char) x;				/* convert from int to char */
			
			if ( (x >= 33 && x < 62) || (x > 62 && x <= 126) ) { 
				*(input + i) = c;	/* place character into memory */
				i += 1;				
			}
			x = fgetc(ifile);	/* read next character */
		}

		if (x == EOF)
			reading = NO;
		printf("Sequence #%d Read. Continuing . . .\n", seq);
		seq += 1;
	}
	printf("\n\nFile Read.\n");

	if (i % NUM_ROW != 0)	 {
		printf("ERROR: SEQEUNCES NOT EQUAL IN LENGTH:\n");
		printf("Proper Memory allocation not possible.");
		printf("CHECK FOR PROPER FASTA STYLE ALIGNMENT WITH GAPS\n");
		exitprogram();
		exit(5);
	}
	else	 {
		NUM_COL = (float) i / NUM_ROW;
		OFFSET = (int) NUM_COL;
	}
	return input;
}





/**************************************************/
/********FUNCTION: score_manager *****************/
/**************************************************/


void	score_manager(char *input)
{
	
	int	**chi_analysis(char *input, int start, int stop);
	
	extern float NUM_COL;
	extern int START;
	extern int STOP;

	int	**chi_matrix;

	int start;	 
	int stop;
	

	start = START;

	if (STOP < 0)	
		stop = (int) NUM_COL;				/*so that program can be run simutaneous on dif. comp*/

	else
		stop = STOP;



	chi_matrix = chi_analysis(input, start, stop);
	
	free_imatrix(chi_matrix, 0, 30, 0, 30);
	 
}





/**************************************************/
/********FUNCTION: chi_analysis *******************/
/**************************************************/

int	**chi_analysis(char *input, int start, int stop)
{

	int		find_unique_elements(char *input, char *aa, int *gap, int var);
	void	crosstab(char *input, int **chi_matrix, char *aa1, char *aa2, int var1, int var2, int count1, int count2);
	
	void 	gap_rectifer(int **chi_matrix,int gap1, int gap2,int *cnt1, int *cnt2);
	void	matrixconverter(int **chi_matrix, int *matrix, double *dmatrix, int count1, int count2, int flag);
	double estExact(int *ROWmatrix, double *dmatrix, double PRECISION, int numrows, int numcols, int *Echeck);
	int		cochrantest(int **chi_matrix, int *rowtot, int *coltot, float *expctd, int numrows, int numcols);
	float	GetFloatTimer (float timer);
	float   chipang(int **nn, int ni, int nj, float *chisq, float *df, double *prob, double *cramrv, float *ccc);

	double *dvector(long nl, long nh);
 
	extern float	NUM_COL;			/*Number of Residue Positions*/
	extern double	PRECISION;			/*precision of calculated probability*/
	extern int		NUM_ROW;			/*Number of Seqeunces*/
	extern int		CELLS;
	extern int		PLACES;



	char	uni1[POSITIONS];			/*array of unique aa in a position*/
	char	uni2[POSITIONS];			/*array of unique aa in a 2nd pos */
	char	*unique1;
	char	*unique2;
	int	count1;							/*number of unique aa in a position*/
	int	count2;							/*number of unique aa in a 2nd pos */
	int	i, j;
	int gap1, gap2;						/*number location of gap (-) in uni arrays*/
	int original_count1;
	/* int test_flag; */
	/** int iqw,jqw,kqw; ***/			/*loop variables for debugging display*/
	
	int *rowtot;						/*array for row totals*/
	int *coltot;						/*array for column totals*/
	float *expctd;						/*matrix for expected frequncy table*/


	float	chi, d, cc;
	double pr,cram;
	float	*chisq, *df, *ccc;
	double  *prob, *cramrv;

	double estprob;						/*probability as determined by exact or estimated exact methods*/
	double logprob;
	int Ecount = 0;
	float Efract = 0.0;
	int *Echeck;						/*Test: has FEXACT ever undergone failure? */
	int Eck = 0;

	float Svalue = 0.0;						/*sensitivity value of scores*/

	int	**chi_matrix;					/*CROSSTABULATION matrix of correlation frequency */
	int **chi_matx2;					/*cross tabulation matrix with origin=ORIGIN */
										/*This ORIGIN shift is necessary for chisquare function*/

	float timer;

	int *matrix_by_rows;				/* 2D chi_matrix represented as 1D array, by rows */
	double *dmatrix;					/* 2D chi_matrix represented as double 1D array, by rows */
	
	char outfile[30] = "alldata";		/*output file for data*/
	char fullname[40];					/*output file for data*/

	char outfile2[30] = "nulldata";		/*output file for meaningless data*/
	char fullname2[40];

	int cramfreq[1001] = {0};
	int index, index2, tt;
	int *lpfreq;

	FILE * output = NULL;				/*output file initiation*/
	FILE * output2 = NULL;
	/* FILE * output3 = NULL;*/


	/**** MEMORY ISSUES ****/


	unique1 = uni1;
	unique2 = uni2;
	chisq = &chi;
	df = &d;
	prob = &pr;
	cramrv = &cram;
	ccc = &cc;
	Echeck = &Eck;

	chi_matrix = imatrix(0, POSITIONS, 0, POSITIONS);     /*allocates memory for an int matrix */ 


	lpfreq = ivector(0,100*prec);	/* allocates memory for log-prob distribution */
 

	for (tt=0; tt <= 100*prec; tt++) {  /*initializes log-prob variable to 0 */
		lpfreq[tt] = 0;
	}



	/*******************************************/
	/*******************************************/
	/*         Open Files for output           */
	/*******************************************/
	/*******************************************/


	sprintf(fullname, "%s%dx%d.txt", outfile, start, stop);
	sprintf(fullname2, "%s%dx%d.txt", outfile2, start, stop);
	

	output = fopen(fullname, "w");
	output2 = fopen(fullname2, "w");


	if (output == NULL) 	 {
		printf("\nFailure to Create alldata File.\n");
		exit(2);
	}

	if (output2 == NULL) 	 {
		printf("\nFailure to Create nulldata File.\n");
		exit(2);
	}


	/*******************************************/
	/***end of opening file for output**********/
	/*******************************************/


	/*******************************************/
	/*******************************************/
	/**********  Start of Main Cycle ***********/
	/*******************************************/
	/*******************************************/

	
		/*******************************************/
		/* Start Timer */
		timer = GetFloatTimer (0.);
		/*******************************************/


	for (i = start; i < stop; i++)	 {

		printf("\n SCORING position #%d's possible interactions. . . please wait.", i);
		
		count1 = find_unique_elements(input, unique1, &gap1, i);					/*finds unique aa in 1st residue*/
		original_count1 = count1;

																				
		for (j = (i+1); j < NUM_COL; j++)  {						
																						 
			/*if (Eck != 0)	{			/*ECheck flag*/
			/*printf("FEXACT has Failed; Monte Carlo Simulation Activated.\n");	*/
			printf("\n   SCORING interaction w/position #%d. Computing.", j);
			/*}*/


			count1 = original_count1; 			/*must be reset each cycle so gap_rectifier does not keep reducing number*/
			count2 = find_unique_elements(input, unique2, &gap2, j);				/*finds unqiue aa in 2nd residue*/

		/**********************/
		/*form crosstab matrix*/
		/**********************/

			crosstab(input, chi_matrix, unique1, unique2, i, j, count1, count2);	/*makes correlation table */
			gap_rectifer(chi_matrix,gap1, gap2, &count1, &count2);					/*adjusts for presence of gaps*/


		/*****************************************************************/
		/* Allocation of memory for reading by various analyitical programs */
		/*****************************************************************/

			rowtot=ivector(0,count1);									/*row totals*/
			coltot=ivector(0,count2);									/*column totals*/
			expctd=vector(0,((count1+1)*(count2+1)-1));					/*expected frequencies*/


		/*****************************************/
		/*Statistical Analysis of crosstab matrix*/
		/*****************************************/

			estprob = MEANINGLESS;	/*meaningless, unless determined */
			*prob = MEANINGLESS;	/*meaningless, unless determined */


			if (count1 >= 1 && count2 >= 1)	{			/*tests that matrix is min. of 2x2*/
				/* test_flag = cochrantest(chi_matrix, rowtot, coltot, expctd, count1+1, count2+1); /*cochran test*/
				/*if (test_flag == 0 || test_flag == 1)	{	/*do normal chi_square test with estimated probability*/
				
				chi_matx2 = subimatrix(chi_matrix,0,POSITIONS,0,POSITIONS,ORIGIN,ORIGIN);/*shifts table to ORIGIN,ORIGIN*/																						  		
				Svalue = chipang(chi_matx2, count1+ORIGIN, count2+ORIGIN, chisq, df, prob, cramrv, ccc);	
																			/*determines chi, cram, df, sensitivity*/	
				estprob = *prob;
				/*}*/

                if (estprob != MEANINGLESS) {        /*do exact or estimated exact probability determination*/
													/*NOTE this second condition is b/c in chisquare, gap percentage check of GAPAMNT see defintions.h*/

					/* allocate memory for different styles of chi matrix */
					matrix_by_rows = ivector(0,((count1+1)*(count2+1)-1));				/*allocate memory*/
					dmatrix = dvector(0,((count1+1)*(count2+1)-1));						/*allocate memory*/
					matrixconverter(chi_matrix, matrix_by_rows, dmatrix, count1, count2, 0);		/*converts 2d matrix to 1d by rows matrix*/
					

					/* call probability calculator */
					estprob = estExact(matrix_by_rows, dmatrix, PRECISION, count1+1, count2+1, Echeck);  /*determine est. exact prob*/
					 
					/* free memory */
					free_subimatrix(chi_matx2, 1, POSITIONS+ORIGIN, 1, POSITIONS+ORIGIN);	/*deallocate memory*/
					free_ivector(matrix_by_rows,0,((count1+1)*(count2+1)-1));					/*deallocate memory*/
				}

				else  { *prob = MEANINGLESS; estprob = MEANINGLESS; } 
			}

			else  {	*prob = MEANINGLESS; estprob = MEANINGLESS;  }		



		/*********************************************************************/
		/* LOG conversion												
		/* counter for sucessful FEXACT function calls						 
		/* FAILURE CHECK for FEXACT analysis							
		/*********************************************************************/


			if ((estprob > 0.0) && (estprob <= 1.0)) 	
						logprob = (-1)*log10(estprob);
			
			else if ((estprob < 0.0) && (estprob != MEANINGLESS))	{	/*EXACT count*/
						estprob = estprob * (-1);
						logprob = (-1)*log10(estprob);
						Ecount = Ecount + 1;
			}
			
			else logprob = 0.0;

		
			
		/* probability using chisquare approximation */

/*			if (pr > 0.0) 	
						logprob = (-1)*log10(pr);

			else logprob = 0.0;

*/

		/*************************/
		/**	OUTPUT to files		**/      /***format: column1, column2, Chi, prob, cram ****/
		/*************************/

				
			if ((pr < 0.0) || (i==j || estprob <0.0) || estprob == 1.0) {
				fprintf(output2, "%d\t%d\t%.3f\t%.5f\t%.3lf\t%.5f\t%f\t%.1f\n", i,j,chi,Svalue,estprob,cram,d,logprob);
			}		/*NOTE: Precision of log prob is 1 place, because calculation of precision is to 10 hits */
					/*SEE function estExact (tally < 10) for further explaination */
			else  {
				
				
					/*************************/
					/* tally distributions.  */
					/*************************/

		
				index = (int) (cram*1000);
				cramfreq[index]++;

				index2 = (int) (logprob*100);
				lpfreq[index2]++;

					/*********************/
				
				
				fprintf(output, "%d\t%d\t%.3f\t%.5f\t%.20lf\t%.5f\t%f\t%.1f\n", i,j,chi,Svalue,estprob,cram,d,logprob);

			}


		/*************************************/
		/**** Free memory for this matrix ****/
		/*************************************/

			free_ivector(coltot,0,(count1));
			free_ivector(rowtot,0,(count2));
			free_vector(expctd,0,((count1+1)*(count2+1)-1));


/*****************************/
/**** Debugging Display ******/
/*****************************/

	/*** start of display ****/

/*
	printf("\n\n\n\nCrosstabulation matrix for %d x %d\n",m,n);
	printf("\nLabels incorrect due to gap rectification\n");

	for (iqw = 0; iqw <= count2; iqw++)
			printf("\t%c",unique2[iqw]);

	printf("\n");

	for (jqw = 0; jqw <= count1; jqw++)	{
			printf("%c",unique1[jqw]);

		for (kqw = 0; kqw <= count2; kqw++)	{
			printf("\t%d", chi_matrix[jqw][kqw]);
		}

		printf("\n");
	}

	printf("\nchi:%f df:%f prob:%f cram:%f cc:%f\n", chi, d, pr, cram, cc);
*/

	/*** END Of DISPLAY  ****/


/***********************************/
/***** END OF DEBUGGING DISPLAY ****/
/***********************************/


 		}
	}
	



	/****************************************/
	/*** END TIMER **************************/
	
	Efract = (float) Ecount;

	timer = GetFloatTimer (timer);
	fprintf(output, "\n\tTotal Elapsed Time is: %.5f\n",timer);
	fprintf(output, "Given: %f EstExact Iterations; %f were EXACT\n",PRECISION, Efract);
	fprintf(output, "Given: %f Columns and %d Rows",NUM_COL, NUM_ROW);
	/****************************************/




	/****close output data pointer file****/
	fclose(output);
	fclose(output2);
	/**************************************/


	return chi_matrix;    /*return pointer for sole purpose of freeing memory allocation*/

}




/**************************************************/
/********FUNCTION: find_unique_elements************/
/**************************************************/


int	find_unique_elements(char *input, char *aa, int *gap, int var)
{

	extern int	OFFSET;
	extern int	NUM_ROW;

	int	i;
	int	j;
	int	k;
	int r;
	int	ext;
	int	counter = 0;		/* # of unique amino acids */


	for (k = 0; k < POSITIONS; k++)	 { 		  /*clear buffer*/
		aa[k] = 0;
	}


	*gap = -1;				/* start off = -1; signifies no gap in position, in any sequence*/

	aa[0] = *(input + var);		/*the first amino acid is always unique*/



	for (i = 1; i < NUM_ROW; i++) {			/*now compare REST of aa=NR-1 :therefore < */
		j = counter;
		ext = NO;

		while (ext == NO) {
			if (*(input + (i * OFFSET + var)) != aa[j]) {

				j = j - 1;
				if (j < 0)		/* No Match, therefore unique */ {
					counter = counter + 1;
					aa[counter] = *(input + (i * OFFSET + var));
					ext = YES;
				}
				else /* No match yet, check previous elements */
				ext = NO;
			}
			else /*MATCHED, so exit*/
			ext = YES;
		}
	}


	for (r = 0; r<=counter; r++) {
		if (aa[r] == '-')			{
			*gap = r;
		}
	}

	return counter;
}


/**************************************************/
/********FUNCTION: Crosstab ***********************/
/**************************************************/

/* Description: using two lists of unique amino acids, found in position var1, var2 */
/* compare all aa in a position to the list [0...p or q]							*/
/* when a match is found in position var, add 1 to matrix position p				*/
/* when a match is found in position var2, add 1 to matrix position q			    */
/* Thus, chi_matrix[p][q] is a tally of frequency of all possible combinations 		*/
/* of amino acids found in var1 with var 2,											*/



void	crosstab(char *input, int **chi_matrix, char *aa1, char *aa2, int var1, int var2, int count1, int count2)
{
	extern int	OFFSET;
	extern int	NUM_ROW;

	void exitprogram();

	int	i, j, k, n, m;
	int p = -1;
	int q = -1;
	int	found_p = 0;
	int	found_q = 0;

	for (j = 0; j <= POSITIONS; j++) {				/*initialize all values to zero */
		for (k = 0; k <= POSITIONS; k++)   {
			chi_matrix[j][k] = 0;
		}
	}

	for (i = 0; i < NUM_ROW; i++)	 {

		n = 0;   		  			/*n are unique aa in 1*/
		/*input compares to n*/
		m = 0;						/*m are unique aa in 2*/
		/*input compares to m*/
		found_p = NO;				/*when a match is found, signal*/
		found_q = NO;				/*when a match is found, signal*/


		while (found_p == NO)   {
			if ( n > count1 )   {
				printf("ERROR: overflow of unique characters set 1 in function crosstab");
				printf("\n n = %d; i = %d",n,i);
				printf("\n %s = aa1 after overflow", aa1);
				printf("\n%c = input",(*input+(i*OFFSET+var1)));
				exitprogram();
				exit(4);
			}
			if (  *(input + (i * OFFSET + var1))  == aa1[n] )   {
				p = n;
				found_p = YES;
			}
			else {
				n += 1;
				found_p = NO;
			}
		}

		while (found_q == NO)   {
			if (m > count2 )	 {
				printf("ERROR: overflow of unique characters set 2 in function crosstab");
				exitprogram();
				exit(3);
			}
			if (  *(input + (i * OFFSET + var2))   == aa2[m] )   {
				q = m;
				found_q = YES;
			}
			else {
				m += 1;
				found_q = NO;
			}
		}

		if (p < 0 || q < 0)	{
			printf("LOGIC ERROR: Character lost in function crosstab");
			exitprogram();
			exit(3);
		}

		chi_matrix[p][q] += 1;

	}
}



/**************************************************/
/**********FUNCTION: gap_rectifier ****************/
/**************************************************/


void gap_rectifer(int **chi_matrix,int gap1, int gap2,int *cnt1, int *cnt2)


/* removes gap (-) from chi_matrix   */


/* NOTE: any zero columns and rows   */
/* produced by gap-as-value removal, */
/* ALSO removed in chisq funct. and EstExact -- REMOVE AS REDUNDANT!!!!!  */

{

	int i,j,k,l,m;
	int rowtot;
	int coltot;

	if (gap1 != -1)  {


		for (i = (gap1+1); i <= *cnt1 ; i++ )   {

			for (j = 0; j <= *cnt2; j++)			{
				chi_matrix[i-1][j] = chi_matrix[i][j];		/*removes gap values*/
			}
		}

		*cnt1 = *cnt1 - 1;

	}


	if (gap2 != -1)	{

		for (k = (gap2+1); k <= *cnt2 ; k++ )   {

			for (l = 0; l <= *cnt1; l++)			{
				chi_matrix[l][k-1] = chi_matrix[l][k];		/*removes gap values*/
			}
		}

		*cnt2 = *cnt2 - 1;

	}



	/* REMOVE ALL ZERO ROWS AND COLUMNS*/


	for (i=0;i <= *cnt1;i++) { 					/*Get the row totals.*/
		rowtot = 0;							/*clears memory*/
		for (j=0;j <= *cnt2;j++) {
			rowtot += chi_matrix[i][j];
		}
		if (rowtot == 0)  {		/*Eliminate any zero rows by reducing the num*/

			for (k = (i+1); k <= *cnt1 ; k++ )   {

				for (m = 0; m <= *cnt2; m++)			{
					chi_matrix[k-1][m] = chi_matrix[k][m];		/*removes gap values*/
				}
			}

			*cnt1 = *cnt1 - 1;
		}


	}


	for (j=0;j <= *cnt2; j++) { 					/*Get the column totals.*/
		coltot=0;							/*clears memory*/

		for (i=0;i<= *cnt1 ;i++) {
			coltot += chi_matrix[i][j];
		}

		if (coltot == 0)   {      /*Eliminate any zero columns.*/

			for (k = (j+1); k <= *cnt2 ; k++ )   {

				for (m = 0; m <= *cnt1; m++)			{
					chi_matrix[m][k-1] = chi_matrix[m][k];		/*removes gap values*/
				}
			}

			*cnt2 = *cnt2 - 1;
		}

	}


}



/**************************************************/
/**********FUNCTION: matrixconverter **************/
/**************************************************/

void matrixconverter(int **chi_matrix, int *matrix, double *dmatrix, int count1, int count2, int flag)
{

	int i,j;
	int index;

	if (flag == 0)	{
	
		for (i = 0; i <= count1; i++)	{
			for (j = 0; j <= count2; j++) {
				index = i*(count2+1) + j;
				matrix[index] = chi_matrix[i][j];
				dmatrix[index] = (double) chi_matrix[i][j];
			}
		}


	}
	else if (flag == 1)	{
	
		for (j = 0; j <= count2; j++)	{
			for (i = 0; i <= count1; i++) {
				index = j*(count1+1) + i;
				matrix[index] = chi_matrix[i][j];
				dmatrix[index] = (double) chi_matrix[i][j];
			}
		}
 
	}

}





/*******************************************/
/**** cochrantest FUNCTION *****************/
/*******************************************/



int cochrantest(int **chi_matrix, int *rowtot, int *coltot, float *expctd, int numrows, int numcols)
{
  
float tally; 
int numrows2,numcols2,i,j,r;
float sum=0.0;


	/***********************/
	/*clears memory = 0    */
	/***********************/


for (r=0;r < numrows*numcols; r++) {
	*(expctd + r) = 0;						/*can probably combine with step below*/
}

	/***********************/
	/*calc row and col info*/
	/***********************/
 
numrows2=numrows; 							/*Number of rows*/
numcols2=numcols; 							/*and columns.*/

for (i=0;i< numrows;i++) { 					/*Get the row totals.*/
	rowtot[i]=0;							/*clears memory*/
	for (j=0;j< numcols;j++) {
		rowtot[i] += chi_matrix[i][j];
		sum += chi_matrix[i][j];
	}
	if (rowtot[i] == 0) --numrows2; 			/*Eliminate any zero rows by reducing the num*/
}


for (j=0;j< numcols;j++) { 					/*Get the column totals.*/
	coltot[j]=0;							/*clears memory*/
	for (i=0;i< numrows;i++) coltot[j] += chi_matrix[i][j]; 
	if (coltot[j] == 0) --numcols2; 			/*Eliminate any zero columns.*/
}


	/**********************************************************/
	/*** Test of Matrix: must be 2x2 after zero removal *******/
	/**********************************************************/

if ((numcols2 < 2) || (numrows2 < 2))	{		/** Chisquare statistical analysis not possible **/
	return MEANINGLESS;
}


	/*****************************************************/
	/* Calculate expected frequencies of original Rmatrix*/
	/*****************************************************/

tally = 0.0;

for (i=0;i< numrows;i++) { 					/*Do the chi-square sum.*/
	for (j=0;j< numcols;j++) {
		*(expctd +(i*numcols + j)) = (coltot[j]*rowtot[i])/sum;

		if (*(expctd +(i*numcols + j)) >= EXPECTED)	
			tally += 1;

		else if (*(expctd +(i*numcols + j)) <= MINIMUM)
			return 1;

	}	
}		


if ( (tally/sum) < PERCENT ) return 1;

else return 0;

 
}
 


/*******************************************/
/**** exitprogram FUNCTION *****************/
/*******************************************/
void exitprogram() {

	float exitprog;
	printf("\n\nPROGRAM EXECUTION COMPLETED");
	printf("\n\n\nPlease enter '1.1' to exit program:");
	scanf("%f",&exitprog);							 
 

	printf("GOOD BYE");

}
/*******************************************/
/**** postionrelater FUNCTION *****************/
/*******************************************/


void positionrelater(char *input)
{

	int i,j;

	
	char	outfile[30] = "Position";		/*output file for data*/
	char	fullname[40];						/*output file for data*/

	extern float NUM_COL;
	extern int OFFSET;
	extern int mainseq;



 
	FILE * output = NULL;						/*output file initiation*/
	sprintf(fullname, "%s.txt", outfile);
	 
	output = fopen(fullname, "w");
 
	if (output == NULL) 	 {
		printf("\nFailure to Create Position File.\n");
		exit(2);
	} 

 
	mainseq = mainseq - 1;	/*because sequences are in data input starting at zero*/
	j = 1;					/*1st character is 1*/


	for(i=0;i < NUM_COL; i++) {

		if ( (*(input + (mainseq * OFFSET + i)) ) != '-' )	{
		
			fprintf(output, "%d\t%d\t%c\n",i,j,(*(input + (mainseq * OFFSET + i))) );
			j = j+1;
		}

	}
 

 
	fclose(output);


}



/*********************************************/
/************  MESSAGE  **********************/
/*********************************************/

void	message(int number)
{

   int  numread;
   int  numwrite;
   char buf[10];
   FILE *inputfile2 = NULL;


   if (number == 1)	 inputfile2 = fopen("scoremsg", "r");/* SCORING MESSAGE */
   if (number == 2)	 inputfile2 = fopen("threshmsg", "r");/* Threshold uses and warnings MESSAGE */
   if (number == 3)	 inputfile2 = fopen("predictmsg", "r");/* Predicition complete MESSAGE */
   if (number == 4)	 inputfile2 = fopen("eliminatemsg", "r");/* eliminate intersecting interactions MESSAGE */
   if (number == 5)	 inputfile2 = fopen("misalignmsg", "r");/* misalignment alg. intro MESSAGE */
   if (number == 6)  inputfile2 = fopen("degenmsg", "r"); /* degenerate predictions MESSAGE */
  
   /* opens input file */
		
	if (inputfile2 == NULL)	printf("\nTEXT NOT FOUND: message %d\n",number);
			   
    else {
		/*read FILE and writes to screen until end of file is reached*/

		while(  !feof(inputfile2)   )
		{
		  numread =  fread(buf,sizeof(char),1,inputfile2); 
		  numwrite = fwrite(buf,sizeof(char),1,stdout);
		} 
   
		fflush(inputfile2);
		fclose(inputfile2);
	}
}

/**************************************************/
/************  FILE_MANAGER  **********************/
/**************************************************/

void	file_manager(void)	{

	/*because progressive threshing -- see Mainscreening.c -- 
	requires a WORK file; it is important that this file is "cleared"
	before the program begins -- otherwise data from previous
	runs will be considered part of the current run.
	For this reason, this function opens then closes the file,
	which has the consequence of clearing it of all data.*/


	char	filename[50] = "PredWORK.txt"; 
	FILE * clearfile = NULL;
	clearfile = fopen(filename, "w");
	fclose(clearfile); 
}



/****************************/
/*!!!! END OF PROGRAM  !!!!!*/
/****************************/

