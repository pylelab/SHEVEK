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

/***********************************************************
The SORTING method below based on code FROM:

Chuck Allison: allison@decus.org or at (801)240-4510.
http://www.freshsources.com/199300f2.htm
modified by Phillip S. Pang: phillip.pang@stanfordalumni.org
************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "shevek.h"
#include "definitions.h"


static int comp(const void *, const void *);
static char tokens [MAXTOKENS] [MAXWIDTH];


/**************************************************************/
/************   redunelim function			*******************/
/**************************************************************/

  
void redunelim(void)	
{
  
	void exitprogram();
	void read_inout(FILE *ifile);
	char	filename[20] = "predWORK.txt";
	char	outfile[25] =  "predictions";		/*output file for data*/
	char	fullname[20];						/*output file for data*/
	char	warning[250];  
	char	star[6] = "666\n";
	char	BLANK[5] = "\0";

	/* variables for sorting *******************/
    size_t i, n;
    size_t idx[MAXTOKENS];
    char *last;
	int tcounter = 1;
	int k;
	/*******************************************/
    size_t strL;

	FILE * inputfile9 = NULL;	 
	FILE * output = NULL;						/*output file initiation*/
 

	/*** open input file ***/
	inputfile9 = fopen(filename, "r");
	if (inputfile9 == NULL)	 {
		printf("\nFAILURE TO OPEN OR LOCATE FILE!!!\n");
		exitprogram();
		exit(99);
	}
	printf("\n\nReading File, Sorting. Please Wait.\n");


	/*** open output file ***/
	sprintf(fullname, "%s.txt", outfile);
	 
	output = fopen(fullname, "w"); 
	if (output == NULL) 	 {
		printf("\nFailure to Create output file.\n");
		exitprogram();
		exit(98);
	}


	
    /*** Read tokens & initialize index array, sorts index ***/
    for (n = 0; fgets(tokens[n],200,inputfile9) != NULL; ++n) {
		printf("\n%s",tokens[n]);
		if ( n >= MAXTOKENS )	{
			printf("\nCRITICAL FAILURE: MEMORY ALLOCATION INSUFFICIENT\n");
			printf("MAXTOKENS definition is too small.\n");
			printf("Contact Programer: phillip.pang@stanfordalumni.org\n");
			exitprogram();
			exit(76);
		}
		else
				idx[n] = n;	
	}

    qsort(idx,n,sizeof idx[0],comp);
	fclose(inputfile9);  /*close input file*/



    /*** Output only unique tokens that are robust/stable ***/

	fprintf(output, "Shevek Predictions\n");
	fprintf(output, "Pos1\tPos2\tCrm V\t-log(P)\tDF\tSensitiv\n"); 


			/* method works because already sorted!! */

	for (k=0;k<250;k++) warning[k] = BLANK[0]; /*flush memory*/ 
 	last = tokens[idx[0]];

    for (i = 1; i < n; ++i)		{
		if ( strcmp(tokens[idx[i]],last) != 0)  {				/*if not equal*/
			if (tcounter >= PSTABILITY )			{					/*if stable*/
				fputs(last, output);
				printf("\n%s",last);
				last = tokens[idx[i]];
				tcounter = 1;
			}
			else	{											/*if not, warn with '*' */			 
				fputs(strcat(strncpy(warning,last,(strlen(last) - 1)),star),output);
				printf("\n**%s",last);
				last = tokens[idx[i]];
				tcounter = 1;
				for (k=0;k<250;k++) warning[k] = BLANK[0]; /*flush memory*/ 
			}						
		}
		else	{ 												/*if equal*/
			tcounter = tcounter + 1;
			printf("\ntcounter = %d", tcounter);
		}
	}	
		   
	last = tokens[idx[n-1]];
	fputs(last, output);

	fclose(output);		/*clost output file*/
      
}

  


/*****************************************************/
/************   COMP function			**************/
/*****************************************************/

 static int comp(const void *p1, const void *p2)	{

    size_t i = * (size_t *) p1;
    size_t j = * (size_t *) p2;

    return strcmp(tokens[i],tokens[j]);
 }


/***************/
/* End of File */
/***************/

