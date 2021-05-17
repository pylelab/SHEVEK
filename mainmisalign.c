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




#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "nrutilp.h"
#include "shevek.h"
#include "definitions.h"

int pfilelength = 0;
int OFFSET2;
int NUM_ROW2;

/********************************************************/
/***** global structures for putative misalignment data */
/********************************************************/

struct position {
	int one;
	int two;
};


struct pair {
	char a;
	char b;
};



struct sequence {
	int name;
	struct position pos;/*position*/
	struct pair cur;	/*current*/
	struct pair sug1;
	struct pair sug2;
	struct pair sug3;
};

struct sequence seq[100000];


/********************************************************/
/********************************************************/
/********************************************************/
	





 
/**************************************************************/
/**************************************************************/
/************   MISALIGN_IDENTIFER function *******************/
/**************************************************************/
/**************************************************************/

/* input is pointer to sequence alignment, offset = OFFSET2, numrows = NUM_ROWS, and pfilelength is
length of prediction.txt file */


void misalign_identifier(char *input, int offset, int numrows)
{
	void	misalignment_manager(char *input,int offset,int numrows);
	void	message(int number);


	extern int OFFSET2;
	extern int NUM_ROW2;
	
	OFFSET2 = offset;
	NUM_ROW2 = numrows;

	message(5);					/*misalignment INTRODUCTORY TEXT */
	
	if (numrows > 100000)	{					/** SEE structure seq declaration **/

		printf("\n\n\n\a\a\a\a\a\aPOSSIBLE MEMORY OVER-RUN in Misalignment Algorithm\n");
		printf(". . .... output is suspect.\n\n\n");
		printf("CONTACT THE PROGRAMER: phillip.pang@stanfordalumni.org");

	}

	
	misalignment_manager(input,offset,numrows);
 	/*Manages multiple matrix files*/
	 

}


/**************************************************************/
/************   MISALIGNMENT_MANAGER function     *************/
/**************************************************************/
   

void	misalignment_manager(char *input,int offset,int numrows)
{

	 
	int	**AR_analysis(char *input, int offset, int numrows, int column1, int column2, int *num_seq);
	int find_unique_seq(int *inp, int *out, int ns);
	double	*openfile4();

	int	**chi_matrix;
	int column1, column2;
	int i,j,tally;

	double	*input4;		/*Predictions Input*/
	int ns = -1;
	int *num_seq;			/*num_seq goes from 0 to num therefore total is actually +1*/
	int num_unique;		/*number of unique sequences*/

	int *inp;				/*1D array of "names"/sequence numbers */
	int *out;				/*array of unique sequences*/ /*see num_unique*/	 


	
	char	outfile[40] = "Misalign";			/*output file for data*/
	char	fullname[100];						/*output file for data*/

	extern struct sequence seq[100000];

	int printflag = NO;

	FILE * output = NULL;						/*output file initiation*/
	num_seq = &ns;

	/** open and read file of prediction data **/

	input4 = openfile4();			/*OPENS THEN READS FILE*/
											/*returns pointer to predictions data*/



	/** initialize output file **/
	sprintf(fullname, "%s.txt", outfile);
	output = fopen(fullname, "w");
 
	if (output == NULL) 	 {
		printf("\nFailure to Create Misalign File.\n");
		exit(2);
	}
 



	/*** FOR EVERY PAIR IN PREDICTION FILE, DO THE FOLLOWING ***/

	for (i=0;i<pfilelength;i++)	{

		/* printf("get this far??"); */
		column1 = (int) (*(input4 + (i * PCOLS + 0)));
		column2 = (int) (*(input4 + (i * PCOLS + 1)));
		/* printf("past first read"); */

		/* Calc. AR values, and identify putative misaligned sequences */
		/* assigns them to structure seq, with seq_num being total number of such sequences */

		chi_matrix = AR_analysis(input, offset, numrows, column1, column2,num_seq);  /*PROBLEM*/

	}



	/*NOW, WITH THIS LIST of putative misaligned sequences, determine which appear frequently*/

	/*printf("past allocation of all seq");*/

	/*** first, convert structure of names to a 1D integer array ***/

	inp = ivector(0,(long) ns);  /*allocate memory; max necessary is num_seq*/
	out = ivector(0,(long) ns);  /*allocate memory; max is num_seq*/   /*IDENTITY OF UNIQUE SEQUENCES*/

	for (i=0;i<=ns;i++)	*(inp+i) = seq[i].name;
	

	
	/*** second, now find all the unique sequences in that array***/
	num_unique = find_unique_seq(inp, out, ns);


	if (num_unique == ns) {			/*both start @ ZERO!!! */		
			printf("\nNO SEQUENCES IDENTIFYABLE AS MISALIGNED\n.");
			fprintf(output,"NO SEQUENCES IDENTIFYABLE AS MISALIGNED\n.");
		}
	
	/*** third, now tally number of times each sequence occurs, output as called for ***/


	else {

		fprintf(output, "SEQUENCES IDENTIFIED AS POSSIBLY MISALIGNED:\n");
		fprintf(output, "(note: sequence numbering starts at 1)\n");
		fprintf(output, "(note: position numbering starts at 0)\n");
		fprintf(output, "(The character '#' denotes NO SUGGESTION POSSIBLE)\n");

		for(i = 0; i <= num_unique; i++)	{
 
			tally = 0;
			for (j= 0; j<= ns; j++)	{				
				if (out[i] == seq[j].name) tally +=1;  
			}


			/* for each unique sequence, if occurs MIN times IS misaligned */
			if (tally >= MINIMUMSEQ)	{


				printflag = YES;

				printf("\n\n\nSequence number %d is possibly misaligned\n", out[i]);
				fprintf(output, "\n\n**********************************************");
				fprintf(output, "\nSequence number %d is possibly misaligned:\n", out[i]);
			


				for (j= 0; j<= ns; j++)	{				
					if (out[i] == seq[j].name)		{  
				
						fprintf(output, "\n\tRegion/Positions:");
						fprintf(output, "%d-%d\t\t", seq[j].pos.one, seq[j].pos.two);
						fprintf(output, "Units/Characters:");
						fprintf(output, "%c-%c\n\n", seq[j].cur.a, seq[j].cur.b);
						fprintf(output, "\tSuggested Alternatives (in order):\n");
						fprintf(output, "\t%c-%c\n", seq[j].sug1.a, seq[j].sug1.b);
						fprintf(output, "\t%c-%c\n", seq[j].sug2.a, seq[j].sug2.b);
						fprintf(output, "\t%c-%c\n", seq[j].sug3.a, seq[j].sug3.b);
					}
				}

			}
		}
	}


	if (printflag == NO) printf("\n\n\nNO SEQUENCES IDENTIFIED");
	free_imatrix(chi_matrix, 0, POSITIONS, 0, POSITIONS);
	fclose(output);

}



/**************************************************************/
/************   AR_ANALYSIS function     **********************/
/**************************************************************/
   

int	**AR_analysis(char *input, int offset, int numrows, int column1, int column2, int *num_seq)
{

	int		find_unique_elements2(char *input, int numrows, int offset, char *aa, int *gap, int var);
	void	crosstab2(char *input, int **chi_matrix, char *aa1, char *aa2, int var1, int var2, int count1, int count2);
	void	arcalc(int **nn, int ni, int nj, float *chisq, float *df, double *prob, float *cramrv, float *ccc, float **ar);

	void 	gap_rectifer1(int **chi_matrix,int gap1, int *cnt1, int *cnt2, char *u1);
	void 	gap_rectifer2(int **chi_matrix,int gap2, int *cnt1, int *cnt2, char *u2);
 
	int cell_finder(char *input,int offset, int numrows, int num_start, char *unique1, char *unique2, int count1, int count2, 
				 int column1, int column2,int **chi_matx2,float **ar_matrix);

	char	outfile[30] = "ARtables";
	char	fullname[40];

	char	uni1[POSITIONS];					/*array of unique aa in a position*/
	char	uni2[POSITIONS];					/*array of unique aa in a 2nd pos */
	char	*unique1;
	char	*unique2;

	int	count1;							/*number of unique aa in a position*/
	int	count2;							/*number of unique aa in a 2nd pos */

	int temp_numseq;
	
	int i,j,k;

	float	chi, d, cram, cc;
	double  pr;
	float	*chisq, *df, *cramrv, *ccc;
	double  *prob;

	int	**chi_matrix;					/*CROSSTABULATION matrix of correlation frequency */
	int **chi_matx2;					/*cross tabulation matrix with origin=ORIGIN */
										/*This ORIGIN shift is necessary for arcalc function*/
										/*can and should be replaced by use of ?vector functions*/

	int gap1, gap2;					/*number location of gap (-) in uni arrays*/


	float **ar_matrix;				/****FOR AR tables!!!! *****/
	
	FILE * output = NULL;

	unique1 = uni1;
	unique2 = uni2;
	chisq = &chi;
	df = &d;
	prob = &pr;
	cramrv = &cram;
	ccc = &cc;

	chi_matrix = imatrix(0, POSITIONS, 0, POSITIONS);     /*allocates memory for an int matrix */

	ar_matrix = matrix(1,POSITIONS+1,1,POSITIONS+1);							/***FOR AR****/

	/*Open File for output, entitled tabout.txt*/

	sprintf(fullname, "%s%dx%d.txt", outfile, column1, column2);

	output = fopen(fullname, "w");

	if (output == NULL) 	 {
		printf("\nFailure to Create File.\n");
		exit(2);
	}



	count1 = find_unique_elements2(input, numrows, offset, unique1, &gap1, column1);
 	count2 = find_unique_elements2(input, numrows, offset, unique2, &gap2, column2);


	/*printf("\nFUE reached for postion %d",column1);*/

	crosstab2(input, chi_matrix, unique1, unique2, column1, column2, count1, count2);

	gap_rectifer1(chi_matrix,gap1, &count1, &count2, unique1);
	gap_rectifer2(chi_matrix,gap2, &count1, &count2, unique2);

	chi_matx2 = subimatrix(chi_matrix,0,POSITIONS,0,POSITIONS,ORIGIN,ORIGIN);

	arcalc(chi_matx2, count1+ORIGIN, count2+ORIGIN, chisq, df, prob, cramrv, ccc, ar_matrix);

	temp_numseq = *num_seq;

	*num_seq = cell_finder(input,offset,numrows,temp_numseq, unique1, unique2, count1, count2, 
				 column1, column2,chi_matx2,ar_matrix);
 
  
	/*************************/
	/*** start of display ****/
	/*************************/

	printf("\n\n\n\nActual Frequency Table for %d x %d\n",column1,column2);

	for (i = 0; i <= count2; i++)
			printf("\t%c",unique2[i]);

	printf("\n");

	for (j = 0; j <= count1; j++)	{
			printf("%c",unique1[j]);

		for (k = 0; k <= count2; k++)	{
			printf("\t%d", chi_matrix[j][k]);
		}

		printf("\n");
	}

	/** printf("\nchi:%f df:%f prob:%.15f cram:%f cc:%f\n", chi, d, pr, cram, cc); **/

	/*************************/
	/*** END Of DISPLAY  ****/
	/*************************/

	/****************************/
	/*** Start of file output ***/
	/****************************/

	fprintf(output,"Actual Frequency Table:");
	fprintf(output, "%d x %d\n", column1,column2);	 


	for (i = 0; i <= (count2); i++)
			fprintf(output, "\t%c",unique2[i]);

	fprintf(output, "\n");

	for (j = 0; j <= (count1); j++)	{
			fprintf(output, "%c",unique1[j]);

		for (k = 0; k <= count2; k++)	{
			fprintf(output, "\t%d", chi_matrix[j][k]);
		}

		fprintf(output, "\n");
	}

	/** fprintf(output, "\nchi:%f df:%f prob:%.15f cram:%f cc:%f\n\n\n", chi, d, pr, cram, cc); **/


	/***output for AR***/
 
	fprintf(output,"\n\nAdjusted Residual (AR) Table:\n");
	for (j = 1; j <= (count1 + 1); j++)	{
		for (k = 1; k <= (count2 + 1); k++)	{
			fprintf(output, "\t%.5f\t", ar_matrix[j][k]);
		}

		fprintf(output, "\n");
	}


	/***END OF OUTPUT FOR AR****/

	fclose(output);

	/*** END of file output***/


	free_subimatrix(chi_matx2, 1, POSITIONS+ORIGIN, 1, POSITIONS+ORIGIN);

	return chi_matrix;    /*return pointer for sole purpose of freeing memory allocation*/

}



/**************************************************************/
/************   FIND_UNIQUE_ELEMENTS2 function     ************/
/**************************************************************/
 


int	find_unique_elements2(char *input, int numrows, int offset, char *aa, int *gap, int var)
{


	int	NUM_ROW2 = numrows;

	int	i;
	int	j;
	int	k;
	int OFFSET2 = offset;	 
	int r;
	int	exit;
	int	counter = 0;		/* # of unique amino acids */


	for (k = 0; k < POSITIONS; k++)	 { 		  /*clear buffer*/
		aa[k] = 0;
	}


	*gap = -1;				/* start off = -1; signifies no gap in position, in any sequence*/

	aa[0] = *(input + var);		/*the first amino acid is always unique*/



	for (i = 1; i < NUM_ROW2; i++) {			/*now compare REST of aa=NR-1 :therefore < */
		j = counter;
		exit = NO;

		while (exit == NO) {
			if (*(input + (i * OFFSET2 + var)) != aa[j]) {

				j = j - 1;
				if (j < 0)		/* No Match, therefore unique */ {
					counter = counter + 1;
					aa[counter] = *(input + (i * OFFSET2 + var));
					exit = YES;
				} else /* No match yet, check previous elements */
					exit = NO;
			} else /*MATCHED, so exit*/
				exit = YES;
		}
	}


	for (r = 0; r<=counter; r++) {
		if (aa[r] == '-')			{
			*gap = r;
			return counter;
		}
	}

	return counter;
}

 
/**************************************************************/
/************   GAP_RECTIFIER1 function     *******************/
/**************************************************************/
 

void gap_rectifer1(int **chi_matrix,int gap1, int *cnt1, int *cnt2, char *u1)


/* removes gap (-) from chi_matrix   */
/* NOTE: any zero columns and rows   */
/* produced by gap-as-value removal, */
/* removed in chisq funct.           */
/* Two gap removers b/c of iterative nature of program */

{

	int i,j;

	if (gap1 != -1)  {


		for (i = (gap1+1); i <= *cnt1 ; i++ )   {

				u1[i-1] = u1[i];							/*removes gap label*/

			for (j = 0; j <= *cnt2; j++)			{
				chi_matrix[i-1][j] = chi_matrix[i][j];		/*removes gap values*/
			}
		}

		*cnt1 = *cnt1 - 1;
	}


}


/**************************************************************/
/************   GAP_RECTIFIER2 function     *******************/
/**************************************************************/


void gap_rectifer2(int **chi_matrix,int gap2, int *cnt1, int *cnt2, char *u2)


{

	int k,l;

	if (gap2 != -1)	{

	for (k = (gap2+1); k <= *cnt2 ; k++ )   {

				u2[k-1] = u2[k];							/*removes gap label*/

			for (l = 0; l <= *cnt1; l++)			{
				chi_matrix[l][k-1] = chi_matrix[l][k];		/*removes gap values*/
			}
		}

		*cnt2 = *cnt2 - 1;
	}

}



/**************************************************************/
/************   OPENFILE4 function     ************************/
/**************************************************************/
 
double	*openfile4()
{

	double	*read_input4(FILE *ifile);
	double	*input4;
	char	filename[50];
 
	FILE * inputfile4 = NULL;

	int	found_file = 0;

	while (found_file == 0) {						/*Queries for input file*/
 
		printf("\nENTER predictions file name (predictions.txt):");

		scanf("%s", filename);
		inputfile4 = fopen(filename, "r");
		if (inputfile4 == NULL)	 {
			printf("\nFAILURE TO OPEN OR LOCATE FILE!!!\n");
			found_file = 0;
		}
		else
		    found_file = 1;
	}

 
 	input4 = read_input4(inputfile4);
	fclose(inputfile4);

	return input4;

}


/**************************************************************/
/************   READ_INPUT4 function     **********************/
/**************************************************************/


double	*read_input4(FILE *ifile)
{
	double *dvector(long nl, long nh);

	void exitprogram();

	double	*input4;

	double inputtemp;

	int i;

	char x;

	extern int pfilelength;



	/*** Memory allocation: determining file size ***/

	printf("\nAllocating Memory for predictions.txt file.\n");


	x = (char) fgetc(ifile);
	while (x != EOF) {

		if (x == '\n')	pfilelength += 1;
 		x = (char) fgetc(ifile);
	}

	pfilelength = pfilelength - 2;  /*see output from threshold.c*/

 
	input4 = dvector(0,pfilelength*PCOLS);
	 
	if (input4 == NULL) {
		printf("\nDynamic Memory Allocation FAILURE!!\n");
		exitprogram();
		exit(1);
	}
	else
		printf("Allocation Completed.\n");
		printf("Allocated memory for %d predictions\n",pfilelength);




	/*** INPUT FILE INTO MEMORY ***/

	fseek(ifile, 0, SEEK_SET);	/*start reading from begining of file*/
	printf("\n\nAttempting to Read Predictions File. Please Wait.\n");

	i = 0; /*start at begining of memory*/
	input4[0] = 0;
	inputtemp = 0;

	fscanf(ifile, "Shevek Predictions\nPos1\tPos2\tCrm V\t-log(P)\tDF\tSensitiv\n");  /*scan through text at top*/
  
	while (i < pfilelength*PCOLS) {
 	
		fscanf(ifile, "%f", &inputtemp);	
		input4[i] = inputtemp;
	
		i = i+1;
	}
 

 	printf("\n\n Predictions File Read.\n");

 	return input4;
}



/**************************************************************/
/************   CELL_FINDER function     **********************/
/**************************************************************/

 
int cell_finder(char *input, int offset, int numrows, int num_start, char *unique1, char *unique2, int count1, int count2, 
				int column1, int column2,int **chi_matx2,float **ar_matrix)	{ 

/*this function returns num_seq --> so that next time it is called, begins there*/	
/*number_start begins at zero, but continues upward*/
/* unique is array from 0 to count */
/* chi_matx2 and ar_matrix go from = to count+1 */
/* column strarts at 0 */

	void sorter(float *A, int left, int right);
	int sequence_finder(char *input, int offset, int numrows, int start, char char1, char char2, int column1, int column2);

	int i;
	int j;
	int temp,temp2,sizetemp;
	int k;
	char charac1;
	char charac2;

	int start;			/*so that if many seq have same character pair, don't return same sequence*/

	int num_seq;		/*counting of number of sequences entered into array*/
	
	float *arvector;

	float best1;
	float best2;
	float best3;

	struct pair s1;
	struct pair s2;
	struct pair s3;

	extern struct sequence seq[100000];
 
	num_seq = num_start;

 /********************************/
 /* main cycle for each AR table */
 /********************************/


	/*************************************************************/
	/* First determine AR values related to misaligned sequences */
	/*************************************************************/

	for (i = 1; i <= (count1+ORIGIN); i++)	{
		for (j = 1; j <= (count2+ORIGIN); j++)	{
	
				/*between 0 and neg threshold*/
			if (   ((ar_matrix[i][j] < 0.0) && (ar_matrix[i][j] >= (ARTHRESH*(-1))) ) 
				|| ((ar_matrix[i][j] > 0.0) && (ar_matrix[i][j] <= ARTHRESH)        )         )	{

				/*is there is a sequence that corrresponds to that cell?*/
				if (chi_matx2[i][j] > 0)	{
					
					temp = chi_matx2[i][j];
					charac1 = unique1[i-1];
					charac2 = unique2[j-1];
					start = 0;

					for (k = 1; k <= temp; k++)	{

						start = sequence_finder(input,offset,numrows,start,charac1,charac2,column1,column2);

						if (start == -1) {
							printf("ERROR ERROR: LOGIC PATHWAY PROBLEM.");
							exit(5);
						}


						num_seq += 1;
						 
						seq[num_seq].name = start;
						seq[num_seq].pos.one = column1;
						seq[num_seq].pos.two = column2;
						seq[num_seq].cur.a = charac1;
						seq[num_seq].cur.b = charac2;

					}
				
				}
					
			}
		
		}
	}



	/********************************************/	
	/* Second, determine suggested alternatives */
	/********************************************/



	
	/*begin by allocating memory for one dimentional array*/
	sizetemp = (count1+ORIGIN)*(count2+ORIGIN);
	arvector = vector(0, sizetemp-1);
	temp2 = 0;

	/*convert two dimentional array into 1 dimentional array*/

	for (i = 1; i <= (count1+ORIGIN); i++)	{
		for (j = 1; j <= (count2+ORIGIN); j++)	{

			arvector[temp2] = ar_matrix[i][j];
			temp2 += 1;

		}
	}

	/*sort array in increasing order*/

	sorter(arvector,0,sizetemp-1);


	/*assign best AR scores to best1...3*/
	/*if less then threshold, then not a best score*/

	if (arvector[sizetemp-1] > ARTHRESH) 	best1 = arvector[sizetemp-1];
	else best1 = (float) NEGONE;

	if (arvector[sizetemp-2] > ARTHRESH) 	best2 = arvector[sizetemp-2];
	else best2 = (float) NEGONE; 
	
	if (arvector[sizetemp-3] > ARTHRESH) 	best3 = arvector[sizetemp-3];
	else best3 = (float) NEGONE;
 

	/*Find characters associated with those scores*/

	for (i = 1; i <= (count1+ORIGIN); i++)	{
		for (j = 1; j <= (count2+ORIGIN); j++)	{

			if (ar_matrix[i][j] == best1)	{
					
					charac1 = unique1[i-1];
					charac2 = unique2[j-1];
					s1.a = charac1;
					s1.b = charac2;

			}

			else if (ar_matrix[i][j] == best2)	{
					
					charac1 = unique1[i-1];
					charac2 = unique2[j-1];
					s2.a = charac1;
					s2.b = charac2;

			}

			else if (ar_matrix[i][j] == best3)	{
					
					charac1 = unique1[i-1];
					charac2 = unique2[j-1];
					s3.a = charac1;
					s3.b = charac2;
			}	
	 
		}
	}

	/* IF as above, best scores were assigned NEGONE, then suggested chacter is '#' */

	if (best1 == NEGONE) {
			s1.a = NOCHAR;
			s1.b = NOCHAR;
	}

	if (best2 == NEGONE) {
			s2.a = NOCHAR;
			s2.b = NOCHAR;
	}

	if (best3 == NEGONE) {
			s3.a = NOCHAR;
			s3.b = NOCHAR;
	}


	/** NOW, assign these suggestions to every sequence from this AR table ***/
	/** obviously, if no identifed sequences, then don't nothing will be assigned**/

	for (k=(num_start+1); k<=num_seq; k++)   {
		
		seq[k].sug1.a =  s1.a;
		seq[k].sug1.b =  s1.b; 
 		seq[k].sug2.a =  s2.a;
		seq[k].sug2.b =  s2.b; 
		seq[k].sug3.a =  s3.a;
		seq[k].sug3.b =  s3.b;

	}



	return num_seq;

}



/**************************************************************/
/************   SEQUENCE_FINDER function     ******************/
/**************************************************************/
 
int sequence_finder(char *input, int offset, int numrows, int start, char char1, char char2, int column1, int column2) {
 
	int i;

	int OFFSET2 = offset;
	int NUM_ROW2 = numrows;
	
	for (i = start; i < NUM_ROW2; i++)	{
		if (  (  (*(input + (i * OFFSET2 + column1))) == char1) && (  (*(input + (i * OFFSET2 + column2))) == char2) )
			return (i+1);
	}
	
	return -1;	/*THIS PATHWAY SHOULD NEVER OCCUR*/
}




/**************************************************************/
/************   FIND_UNIQUE_SEQ function   ********************/
/**************************************************************/
 
int	find_unique_seq(int *inp, int *out, int num_seq)
{

	int	i;
	int	j;
	int	k;
 
	int	exit;
	int	counter = 0;		/* # of unique sequences */

 

	for (k = 0; k <= num_seq; k++)	  		  /*clear buffer*/
		out[k] = 0;								/* need this??*/
	

	out[0] = *(inp + 0);		/*the first sequence is always unique*/


	for (i = 1; i <= num_seq; i++) {			/*now compare REST */
		j = counter;
		exit = NO;

		while (exit == NO) {
			if ( *(inp + i) != out[j]) {

				j = j - 1;
				if (j < 0)		/* No Match, therefore unique */ {
					counter = counter + 1;
					out[counter] = *(inp + i);
					exit = YES;
				} else /* No match yet, check previous elements */
					exit = NO;
			} else /*MATCHED, so exit*/
				exit = YES;
		}
	}

	return counter;
}

/**************************************************/
/********FUNCTION: Crosstab2 ***********************/
/**************************************************/

/* Description: using two lists of unique amino acids, found in position var1, var2 */
/* compare all aa in a position to the list [0...p or q]							*/
/* when a match is found in position var, add 1 to matrix position p				*/
/* when a match is found in position var2, add 1 to matrix position q			    */
/* Thus, chi_matrix[p][q] is a tally of frequency of all possible combinations 		*/
/* of amino acids found in var1 with var 2,											*/

/** THIS FUNCTION IS REPEATED because it calls two externals, that are limited to another file**/


void	crosstab2(char *input, int **chi_matrix, char *aa1, char *aa2, int var1, int var2, int count1, int count2)
{
	extern int	OFFSET2;
	extern int	NUM_ROW2;

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

	for (i = 0; i < NUM_ROW2; i++)	 {

		n = 0;   		  			/*n are unique aa in 1*/
		/*input compares to n*/
		m = 0;						/*m are unique aa in 2*/
		/*input compares to m*/
		found_p = NO;				/*when a match is found, signal*/
		found_q = NO;				/*when a match is found, signal*/


		while (found_p == NO)   {
			if ( n > count1 )   {
				printf("ERROR: overflow of unique characters set 1");
				printf("\n n = %d; i = %d",n,i);
				printf("\n %s = aa1 after overflow", aa1);
				printf("\n%c = input",(*input+(i*OFFSET2+var1)));
				exitprogram();
				exit(4);
			}
			if (  *(input + (i * OFFSET2 + var1))  == aa1[n] )   {
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
				printf("ERROR: overflow of unique amino acids set 2");
				exitprogram();
				exit(3);
			}
			if (  *(input + (i * OFFSET2 + var2))   == aa2[m] )   {
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







/*********************************************/
/*********************************************/
/*********************************************/
 
/* standlib functions as defined from ANSI C */

/*********************************************/
/*********************************************/
/*********************************************/

  

/**************************************************************/
/************   SORTER function   *****************************/
/**************************************************************/

void sorter(float *A, int left, int right) {

	int i, last;
	void swapper (float *A, int i, int j);
	void sorter(float *A, int left, int right);  /*redundant prototype for recursion*/
	
	if (left >= right)
		return;
	swapper(A, left, (left+right)/2);
	last = left;
	
	for (i=left+1; i <=right; i++)	
		if ( (*(A+i)) < (*(A+left))  )    swapper(A, ++last, i);
		 

	swapper(A,left,last);
	sorter(A,left, last-1);
	sorter(A,last+1,right);
}


/**************************************************************/
/************   SWAPPER function   ****************************/
/**************************************************************/



void swapper(float *A, int i, int j)	{

	float temp;

	temp = (*(A+i));
	(*(A+i)) = (*(A+j));
	(*(A+j)) = temp;

}

