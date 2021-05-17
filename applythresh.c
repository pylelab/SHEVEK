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
 
#include "nrutilp.h" 
#include "shevek.h"
#include "definitions.h"


double POFF;
double CRV_CUTOFF;
double SOFF;
double DFOFF;
float RES;

int FILEL1;
int FILEL2;
  

/**************************************************************/
/************   APPLY_THRESHOLDS function   *******************/
/**************************************************************/

  
void apply_thresholds(void)
{
 
	double	*openfile3(int pass,float *actpmax);
	int		prelimscan(double *input, int filelength); 
	void	eliminator2(double *input,int filelength);
	void	exitprogram();
	void	free_dvector(double *v, long nl, long nh);
	void	redunelim(void);
	void	message(int number);

	double	*input;
	double *input2;
 
	float *actpmax;
	float apm = 0.0;
 
	int filelength2;

	int file1max = 0;
	int file2max = 0;

	float sdef;
	float dfdef;
	float resdef;

	actpmax = &apm;

	/* OPEN FILE for input reading*/
 	input = openfile3(0,actpmax);							/*OPENS THEN READS FILE*/
	
	printf("actpmax = %f",apm);

	/*returns pointer to sequence data*/
	if (input[0] == MEANINGLESS) {
		printf("INPUT DATA FILE IS EMPTY -- Error!! exiting . . .");
		exitprogram();
		exit(10);

	}
  


	/***************** obtains parameters from user ************************/

	sdef = (float) SDEFAULT;
	dfdef = (float) DFMAX;
	resdef = (float) RESDEFAULT;
	

	printf("\n***********************************************\n");
	printf("PARAMETER ENTRY (use numbers above as defaults)\n");
	printf("***********************************************\n");

	/*printf("\nEnter -log(P) resolution (.1 to 1;default:%.2f): ", resdef);		
	scanf("%f",&RES);	*/
	RES = resdef;

	printf("\nPlease enter a -log(P) threshold (<=: range:0-inf): ");		
	scanf("%f",&POFF);								/*query for p-value parameter*/
													 

	printf("\nPlease enter the default V threshold (<=: range:0-1): ");
	scanf("%f",&CRV_CUTOFF);						/*query for Cramer's V paramerter*/
												 

	printf("\nPlease enter a Sensitivity threshold (<: default:%.1f): ", sdef); 
	scanf("%f",&SOFF);		
	
	
	printf("\nPlease enter a df threshold (<: default:%.0f): ", dfdef); 
	scanf("%f",&DFOFF);		
		
 


	/************* DETERMINE Predictions of progressive logP CUTOFFS ***************/
	message(4);	
	while (POFF < apm)	{
 	
		filelength2 = prelimscan(input,FILEL1);			/*Outputs File PVSDscan.txt */ 									 
		input2 = openfile3(1,actpmax);					/*Opens file just created -- PVSDscan.txt*/

		if (input2[0] != MEANINGLESS)	{				/*see read_input3() -- file is not empty*/
			eliminator2(input2,FILEL2);					/*ELIMINATE INTERSECTING INTERACTIONS*/											
			if (FILEL1 > file1max) file1max = FILEL1;   /* get maxfile length to free memory */
			if (FILEL2 > file2max) file2max = FILEL2;
		}
		POFF = POFF + RES;
	}

	POFF = apm;										/* final loop for highest -logP value */
	filelength2 = prelimscan(input,FILEL1);			/*Outputs File PVSDscan.txt */ 									 
	input2 = openfile3(1,actpmax);					/*Opens file just created -- PVSDscan.txt*/

	if (input2[0] != MEANINGLESS)	{				/*see read_input3() -- file is not empty*/
		eliminator2(input2,FILEL2);					/*ELIMINATE INTERSECTING INTERACTIONS*/											
		if (FILEL1 > file1max) file1max = FILEL1;   /* get maxfile length to free memory */
		if (FILEL2 > file2max) file2max = FILEL2;
	}

	/*****************************************************************************/


	/************** FREE MEMORY ASSOCIATED WITH PROGRESSIVE SCANS ****************/ 
	free_dvector(input,0,file1max*COLS);	/*frees memory allocated for input file	*/											
	free_dvector(input2,0,file2max*COLS);





	/************** ELIMINATE FROM FILE "predWORK" any REDUNDANT PREDICTIONS****************/
	redunelim();

}




/**************************************************************/
/************   OPENFILE3 function   **************************/
/**************************************************************/
 
/* NOTE openfile2 was combined here; no longer exists as function*/

double	*openfile3(int pass, float *actpmax)
{

	double	*read_input3(FILE *ifile, int pass, float *actpmax);
	double	*input;
	char	filename[50];

	FILE * inputfile = NULL;

	int	found_file = 0;

	while (found_file == 0) {						/*Queries for input file*/

		if (pass == 0)	{

		printf("\n\nEnter alldata file name for theshing (alldataXX.txt):");
		scanf("%s", filename);
														 
	}

 
		if (pass == 1) strcpy(filename, "PVSDscan.txt");
 
		
		inputfile = fopen(filename, "r");
		if (inputfile == NULL)	 {
			printf("\nFAILURE TO OPEN OR LOCATE FILE!!!\n");
			found_file = 0;
		}
		else
		    found_file = 1;
	}

 
 	input = read_input3(inputfile, pass, actpmax);
	
	fclose(inputfile);

	return input;

}



/**************************************************************/
/************   READ_INPUT3 function   ************************/
/**************************************************************/
 
/*note read_input2 was combined here; no longer exists as function*/

double	*read_input3(FILE *ifile, int pass, float *actpmax)
{
	void exitprogram();

	double *dvector(long nl, long nh);
	double	*input;

	double inputtemp;
	int i;
	int x;
	int filelength =0;
 
/*** Memory allocation: determining file size ***/

	printf("\nAllocating Memory for data file.\n");
 
	x = fgetc(ifile);
	while (x != EOF) {

		if (x == '\n')	filelength += 1;
 		x = fgetc(ifile);
	}

	if (filelength > 0)	{

		if (pass == 0)	{	/* alldata file */
			filelength = filelength - 3;  /** see chi_analysis for output file extras of alldata.txt **/

			FILEL1 = filelength;
		}

		if (pass == 1 || pass == 2) {	/*PVSDscan file or allpredict.txt */

			FILEL2 = filelength;
		}


		
		input = dvector(0,filelength*COLS);

		 
		if (input == NULL) {
			printf("\nDynamic Memory Allocation FAILURE!!\n");
			exitprogram();
			exit(1);
		}
		else
			printf("Allocation Completed.\n");


 
		/*** INPUT FILE INTO MEMORY ***/

		fseek(ifile, 0, SEEK_SET);	/*start reading from begining of file*/
		printf("Attempting to Read File. Please Wait.\n");

		i = 0; /*start at begining of memory*/
		input[0] = 0;
		inputtemp = 0;


 
		while (i < filelength*COLS) {
 		
			fscanf(ifile, "%f", &inputtemp);	
			input[i] = inputtemp;

			/* find best -logP */
			if (    (pass == 0) 
				&&  ( i % COLS == col2 )
				&&  (input[i] > *actpmax)      )				{ /*i.e if in -logP column*/
				*actpmax = (float) input[i];
			}
			i = i+1;
		}
 

 		printf("\n\nFile Read.\n");
	}

	else  {
		printf("NOTE: EMPTY FILE! No Predictions for this threshold: %.1lf\n", POFF);
		input = dvector(0,10);  /*RANDOM DUMB ALLOCATION OF MEMORY*/
		input[0] = MEANINGLESS;
	}

 	return input;
}


/**************************************************************/
/************   PRELIMSCAN function     ***********************/
/**************************************************************/


int	prelimscan(double *input, int filelength)
{


extern double POFF;
extern double CRV_CUTOFF;	
extern double SOFF;
extern double DFOFF;

int filelength2 = 0;



char	outfile[40] = "PVSDscan";		/*output file for data*/
char	fullname[50];						/*output file for data*/

 
FILE * output = NULL;						/*output file initiation*/
 

/*i = row and col = col*/
int i=0;


sprintf(fullname, "%s.txt", outfile);
	 

output = fopen(fullname, "w");
 
if (output == NULL) 	 {
	printf("\nFailure to Create PVSDscan File.\n");
	exit(2);
}
 
 
for (i = 0; i < filelength; i++)	{
	
	/* If thresholds are met, include in PVSDscan file */

	if (      ( (*(input + (i * COLS + col1) ) ) >= CRV_CUTOFF)			/*see definitions.h*/
		  &&  ( (*(input + (i * COLS + col2))) >= POFF)          
		  &&  ( (*(input + (i * COLS + col5))) < SOFF)       
		  &&  ( (*(input + (i * COLS + col3))) < DFOFF)			) {
 
			filelength2 = filelength2 + 1;

			fprintf(output, "%d\t%d\t%.3f\t%.5f\t%.20lf\t%.5f\t%f\t%.5f\n", 
				(int) (*(input + (i * COLS + 0))),
				(int) (*(input + (i * COLS + 1))),
				(float) (*(input + (i * COLS + 2))),
				(float) (*(input + (i * COLS + 3))),
				(*(input + (i * COLS + 4))),
				(float) (*(input + (i * COLS + 5))),
				(float) (*(input + (i * COLS + 6))),
				(float) (*(input + (i * COLS + 7)))  ); 

	}

 
}
 
 
fclose(output);

return filelength2;

}


