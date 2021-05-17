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



#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
 
#include "nrutilp.h" 
#include "shevek.h"
#include "definitions.h"

double *dinput; 

/**************************************************************/
/************   DISTRIBUTION_ANALYZER function       **********/
/**************************************************************/

	static int pcomp(const void *p1, const void *p2);
	static int vcomp(const void *p1, const void *p2);
	static int dfcomp(const void *p1, const void *p2);

void distribution_analyzer(float *VThresh, float *PThresh, float *DFThresh)
{ 
	double *dvector(long nl, long nh);
	int *ivector(long nl, long nh);
	void exitprogram();
    /*
	static int pcomp(const void *p1, const void *p2);
	static int vcomp(const void *p1, const void *p2);
	static int dfcomp(const void *p1, const void *p2);
*/

	unsigned int i, n;
	unsigned int filelength;

	int x;

	double inputtemp;

	int		*pidx;
	int		*vidx;
	int		*dfidx;
	
    double	Viqr, Piqr, DFiqr;
	double	pmedian, vmedian, dfmedian;
	double	firstqV,thirdqV,firstqP,thirdqP,firstqDF, thirdqDF;

	double	 BC;  /*Bonforroni Correction*/
	 
	
	char	filename[50];

	FILE * inputfile = NULL;
	int	found_file = 0;
 
 
/***************** OPEN FILE ******/

	while (found_file == 0) {						/*Queries for input file*/
		printf("\nENTER alldata file name for distribution analysis:");
		scanf("%s", filename);
		inputfile = fopen(filename, "r");
		if (inputfile == NULL)	 {
			printf("\nFAILURE TO OPEN OR LOCATE FILE!!!\n");
			found_file = 0;
		}
		else
		    found_file = 1;
	}



/**************** READ FILE ********/

	printf("\nAllocating Memory for file.\n"); 
	
	/** determine file size **/		

	filelength = 0;
	x = fgetc(inputfile);
	
	while (x != EOF) {
		if (x == '\n')	filelength += 1;
 		x = fgetc(inputfile);
	}
	filelength = filelength - (int) ALLOFFSET;


	/*Allocate memory */	
	dinput = dvector(0,filelength*COLS);
	 
	if (dinput == NULL) {
		printf("\nDynamic Memory Allocation FAILURE!!\n");
		exitprogram();
		exit(1);
	}
	else
		printf("Allocation Completed.\n");

 
	/*** INPUT FILE INTO MEMORY ***/
	fseek(inputfile, 0, SEEK_SET);	/*start reading from begining of file*/
	printf("Attempting to Read File. Please Wait.\n");
	
	i = 0; /*start at begining of memory*/
	dinput[0] = 0;
	inputtemp = 0;

	while (i < filelength*COLS) {		
			fscanf(inputfile, "%f", &inputtemp);	
			dinput[i] = inputtemp;
			i = i+1;
	}
 
 	printf("\n\nAlldata File Read.\n");
	fclose(inputfile);





/****************** SORTING -logP, V, df DATA *****/

	/* allocate memory for indexes */

	pidx = ivector(0, filelength-1);
	vidx = ivector(0, filelength-1);
	dfidx = ivector(0, filelength-1);
	

	/*initialize indexes */
	for (n = 0; n < filelength; ++n) {
		pidx[n]  = n;	
		vidx[n]  = n;
		dfidx[n] = n;
	}

	/*sorts indexes */

	printf("\nSORTING -logP,V, df values for analysis . . .\n");
    qsort(pidx,n,sizeof pidx[0],pcomp);		/* sort -logP index */
	qsort(vidx,n,sizeof vidx[0],vcomp);		/* sort V index */
	qsort(dfidx,n,sizeof dfidx[0],dfcomp);		/* sort df index */




/****************** FIND MEDIAN AND QUARTILE VALUES *****/

   	/* median */
	if (filelength % 2 != 0)	{
			pmedian = (   (*(dinput + ((pidx[filelength/2]) * COLS + col2)))  
				        + (*(dinput + ((pidx[(filelength/2) + 1]) * COLS + col2)))  )/2;
			vmedian = (   (*(dinput + ((vidx[filelength/2]) * COLS + col1)))  
				        + (*(dinput + ((vidx[(filelength/2) + 1]) * COLS + col1)))  )/2;
			dfmedian = (  (*(dinput + ((dfidx[filelength/2]) * COLS + col3)))  
				        + (*(dinput + ((dfidx[(filelength/2) + 1]) * COLS + col3)))  )/2;
	
	
	}
	else	{
			pmedian = (*(dinput + ((pidx[filelength/2]) * COLS + col2)));
			vmedian = (*(dinput + ((vidx[filelength/2]) * COLS + col1)));
			dfmedian = (*(dinput + ((dfidx[filelength/2]) * COLS + col3)));
	}


	/* 1st quartile */
	if (filelength % 4 != 0)	{
 			firstqP = (   (*(dinput + ((pidx[filelength/4]) * COLS + col2)))  
				        + (*(dinput + ((pidx[(filelength/4) + 1]) * COLS + col2)))  )/2;
			firstqV = (   (*(dinput + ((vidx[filelength/4]) * COLS + col1)))  
				        + (*(dinput + ((vidx[(filelength/4) + 1]) * COLS + col1)))  )/2;
			firstqDF = (  (*(dinput + ((dfidx[filelength/4]) * COLS + col3)))  
				        + (*(dinput + ((dfidx[(filelength/4) + 1]) * COLS + col3)))  )/2;
	}
	else	{
			firstqP = (*(dinput + ((pidx[filelength/4]) * COLS + col2)));
			firstqV = (*(dinput + ((vidx[filelength/4]) * COLS + col1)));
			firstqDF = (*(dinput + ((dfidx[filelength/4]) * COLS + col3)));
	}


	/* third quartile*/
	if ((3*filelength) % 4 != 0)	{
 			thirdqP = (   (*(dinput + ((pidx[3*filelength/4]) * COLS + col2)))  
				        + (*(dinput + ((pidx[(3*filelength/4) + 1]) * COLS + col2)))  )/2;
			thirdqV = (   (*(dinput + ((vidx[3*filelength/4]) * COLS + col1)))  
				        + (*(dinput + ((vidx[(3*filelength/4) + 1]) * COLS + col1)))  )/2;
			thirdqDF = (  (*(dinput + ((dfidx[3*filelength/4]) * COLS + col3)))  
				        + (*(dinput + ((dfidx[(3*filelength/4) + 1]) * COLS + col3)))  )/2;
	}
	else	{
			thirdqP = (*(dinput + ((pidx[3*filelength/4]) * COLS + col2)));
			thirdqV = (*(dinput + ((vidx[3*filelength/4]) * COLS + col1)));
			thirdqDF = (*(dinput + ((dfidx[3*filelength/4]) * COLS + col3)));
	}


	/* calculate INTERQUARTILE RANGE */
	Piqr = thirdqP - firstqP;
	Viqr = thirdqV - firstqV;
	DFiqr = thirdqDF - firstqDF;

	 
	/* Calculation of Bonforroni Correction */
	BC = (0.01)/((float) filelength);
	BC = (-1)*log(BC);

	
	/* output to screen calculated values */
     
	printf("\n-log(P) Median = %.1f\n",pmedian);
	printf("-log(P) IQR = %.1f\n",Piqr);
	printf("-log(P) 1Q-value = %.1f\n",firstqP);
	printf("-log(P) 3Q-value = %.1f\n",thirdqP);
	printf("-log(p) Bonforroni = %.1f\n",BC); 

	printf("\nV Median = %.3f\n",vmedian);
	printf("V IQR = %.3f\n",Viqr);
	printf("V 1Q-value = %.3f\n",firstqV);
	printf("V 3Q-value = %.3f\n",thirdqV);

	printf("\ndf Median = %.1f\n",dfmedian);
	printf("df IQR = %.1f\n",DFiqr);
	printf("df 1Q-value = %.1f\n",firstqDF);
	printf("df 3Q-value = %.1f\n",thirdqDF);


	/*Calculate Thresholds*/
 
	*VThresh = (float) (thirdqV + (VMULT * Viqr));
	*PThresh = (float) (thirdqP + (PMULT * Piqr));
	*DFThresh = (float) firstqDF;


	/*******************************************/
	/* Threshold Warnings/Triggers/Adjustments */


	if (*VThresh < (float) VMINIMUM)		{
			*VThresh = (float) VMINIMUM;
			printf("\n\tWARNING: Vminimum Triggered!\n");
			printf("\tDataset may contain too little information\n");
	}

	
	if ( (*VThresh >= (float) VMINIMUM) && (*VThresh <= (float) VMAXIMUM) )		{
			*VThresh = (float) VMINIMUM;		/*V thresh should not be used for elimination*/
	}	



	if (*VThresh > (float) VMAXIMUM)		{
			printf("\n\tWARNING: Vmaximum Triggered!\n");
			printf("\tDataset may contain too little diversity\n");
			printf("\tNo predictions possible at this threshold.\n");
			printf("\t Hit CNTR-C to end program.\n");
	}			


	if (*PThresh < (float) PMINIMUM)		{
			*PThresh = (float) PMINIMUM;
			printf("\n\tWARNING: Pminimum Triggered!\n");
			printf("\tDataset may contain too few sequences\n");
	}



	/************************/

}

 


/*****************************************************/
/************   pCOMP function			**************/
/*****************************************************/

 static int pcomp(const void *p1, const void *p2)	{

    size_t i = * (size_t *) p1;
    size_t j = * (size_t *) p2;

	if		(  (*(dinput + (i * COLS + col2)))  <  (*(dinput + (j * COLS + col2)))   )	return -1;
	else if (  (*(dinput + (i * COLS + col2))) == (*(dinput + (j * COLS + col2)))   )		return 0;
	else	return 1;
 }




/*****************************************************/
/************   vCOMP function			**************/
/*****************************************************/

 static int vcomp(const void *p1, const void *p2)	{

    size_t i = * (size_t *) p1;
    size_t j = * (size_t *) p2;
	
	if		(  (*(dinput + (i * COLS + col1)))  <  (*(dinput + (j * COLS + col1)))   )	return -1;
	else if (  (*(dinput + (i * COLS + col1))) == (*(dinput + (j * COLS + col1)))   )		return 0;
	else	return 1;
 }
 

 
/*****************************************************/
/************   dfCOMP function			**************/
/*****************************************************/

 static int dfcomp(const void *p1, const void *p2)	{

    size_t i = * (size_t *) p1;
    size_t j = * (size_t *) p2;
	
	if		(  (*(dinput + (i * COLS + col3)))  <  (*(dinput + (j * COLS + col3)))   )	return -1;
	else if (  (*(dinput + (i * COLS + col3))) == (*(dinput + (j * COLS + col3)))   )		return 0;
	else	return 1;
 }

/*************/
/*END OF FILE*/
/*************/
