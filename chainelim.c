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



/*************************************************/
/************  ELIMINATOR2  **********************/
/*************************************************/
/*historical note:								 */
/*eliminator 1 was an older method of elimination*/
 

void	eliminator2(double *input,int filelength)
{

	void	pathfinder(int fileposition, int group, double *datafile, int filemax);
	void	message(int number);
	void	logicerror(int code);
	void	exitprogram();

  
	int		group = -1;

	float	dferror;
	float	dftemp;
	
	int		groupticker;
	int		groupnumber;

	int		i;

	double	crambest = 0;
	int		firstv;
	
	int		firstp;
	double	lprobbest = 0;

	int		firstdf;
	int		dfbest = 0;

	int		printcounter = 0;

	char	outfile[40] = "PredWORK";		/*output file for data*/
	char	fullname[50];						/*output file for data*/

	FILE * output = NULL;						/*output file initiation*/
   





	/* Part One -- FIND ALL NETWORKS of INTERACTRIONS/REDUNDANT INTERACTIONS*/
	/*for every postion not found already as part of a network, find its network, and assign own group name */


	for (i=0;i<filelength;i++) {

		if (  (*(input + (i * COLS + col4))) >= 0    )	{	/*not assigned to a group*/

			pathfinder(i, group, input, filelength);				/*find its network, and assign a group name to entire network*/

			group = group - 1;											/*get new group name*/

		}
	}







	/* Part Two -- Having assigned a group to every preliminary association interaction, find bestin each group */
	/* best V to error of .05; then best P; then lowest DF */


	group = group + 1;		/*since last assignment will decrement */
	
	groupticker = group;


	/* for debugging */

	/* is this bug source??  *

/*	if (groupticker = -1)	{
		groupticker = groupticker;
	}
*/ 


	while (groupticker < 0)	{			/*for each and every group, find best V */

		firstp = NO;
		groupnumber = 0;


/*find best P*/
/*NOTE for P, BEST is HIGHEST*/


		for (i=0;i<filelength;i++) {
			

			if ( (*(input + ( i * COLS + col4))) == groupticker ) {

				groupnumber = groupnumber + 1;  /*count number of members in the group*/
				
			}


			if ( firstp == NO && (*(input + ( i * COLS + col4))) == groupticker  )	{   /*find first member of this group */

				lprobbest = (*(input + (i * COLS + col2)));		/*assign best P to first P in group*/
				firstp = YES;

			}

			else if ( firstp == YES && (*(input + ( i * COLS + col4))) == groupticker  )	{   /*find other members of this group */

				if ( lprobbest  < (*(input + ( i * COLS + col2)))     )		{       /*compare P's within group*/   
						lprobbest = (*(input + ( i * COLS + col2)));		/*if new one is higher, reasign best */									
					}
			}

		}

		/*eliminate all but best P, within Perror*/

		if (lprobbest == 0) {
			printf("LOGIC ERROR: no best P in file, function eliminator2");
			exitprogram();
			exit(5);
		}

		for (i=0;i<filelength;i++) {		
			
			if (     (*(input + ( i * COLS + col4))) == groupticker							/*same group*/
				&&   (*(input + ( i * COLS + col2)))  <  (lprobbest - ((float) PERROR))	)	{		/*less then best + error*/

				*(input + (i * COLS + col1)) = NEGONE;	
				groupnumber = groupnumber - 1;

			}
		}

/*is there more than one best left?*/
/*then judge by df */
/*NOTE for df, BEST is lowest */


		if (groupnumber > 1) {			/*i.e. more than one member in the group left -- need to use logP as judge */

			firstdf = NO;

/*find best DF*/

			for (i=0;i<filelength;i++) {		
				
				if (     (*(input + ( i * COLS + col4))) == groupticker						/*find first member of group and sub group*/
					&&   (*(input + ( i * COLS + col1)))  != NEGONE							/*sub-group of best V, within error*/
					&&	  firstdf == NO									) {					 


					dfbest = (int) (*(input + (i * COLS + col3)));		/*assign best df to first df in group and sub group*/
					firstdf = YES;

				}


				else if (       firstdf == YES										/*find other members of this group */
					       &&   (*(input + ( i * COLS + col4))) == groupticker  
						   &&	(*(input + ( i * COLS + col1)))  != NEGONE	) {   


						if ( dfbest  > ((int) (*(input + ( i * COLS + col3))))     )		{       /*compare DF's within group*/   
							dfbest = (int) (*(input + ( i * COLS + col3)));		/*if new one is lower, reasign best */									
						}
				}



			}

/*eliminate all but best*/

	
			if (dfbest == 0) {
				printf("LOGIC ERROR: no best df in file, function eliminator2");
				exitprogram();
				exit(5);
			}

			/* calculation of DF error */

			dftemp = (float) sqrt( ((double)dfbest) );
			dftemp = dftemp + ((int) DFDIMERR);
			dferror = dftemp*dftemp;

			/***************************/


			for (i=0;i<filelength;i++) {		
				
				if (      (*(input + ( i * COLS + col4))) == groupticker							/*same group*/
					 &&	  (*(input + ( i * COLS + col1)))  != NEGONE							/*and not already elim*/
					 &&   (*(input + ( i * COLS + col3)))  >  (dfbest + dferror)	)	{					/*greater then df best*/

						*(input + (i * COLS + col1)) = NEGONE;	
						groupnumber = groupnumber - 1;

				}
			}

		}


		if (groupnumber > 1)	{		/*i.e. after DF and P best, still more than one, judge by V*/

/*is there more than one member in the group left?*/
/*then judge by df */

			firstv = NO;


/*find best df and eliminatring rest as go*/
/*NOTE for V, BEST is HIGHEST*/

			for (i=0;i<filelength;i++) {		
				
				if (     (*(input + ( i * COLS + col4))) == groupticker						/*find first member of group and sub group*/
					&&   (*(input + ( i * COLS + col1)))  != NEGONE							/*sub-group */
					&&	  firstv == NO									) {					 


					crambest = (*(input + (i * COLS + col1)));		/*assign best to first in group and sub group*/
					firstv = YES;

				}


				else if (       firstv == YES										/*find other members of this group */
					       &&   (*(input + ( i * COLS + col4))) == groupticker  
						   &&	(*(input + ( i * COLS + col1)))  != NEGONE	) {   


						if ( crambest  <  (*(input + ( i * COLS + col1)))     )		{       /*compare within group*/   
							crambest = (*(input + ( i * COLS + col1)));		/*if new one is lower, reasign dfbest */									
						}
				}

			}	/*loop for finding best V by line*/
 


			/*Next, eliminate all but best V, within error*/

			if (crambest == 0) {
				printf("LOGIC ERROR: no best V in file, function eliminator2");
				exitprogram();
				exit(5);
			}


			for (i=0;i<filelength;i++) {		
				
				if (     (*(input + ( i * COLS + col4))) == groupticker							/*same group*/
					&&	 (*(input + ( i * COLS + col1)))  != NEGONE								/*and not already elimin*/
					&&   (*(input + ( i * COLS + col1)))  <  crambest	)	{		/*less then best + error*/

					*(input + (i * COLS + col1)) = NEGONE;	
					groupnumber = groupnumber - 1;

				}
			}


			/*Next, determine if group still has more than one member; if so indicate as NEGATIVE*/

			if (groupnumber > 1)	{
	
				message(6);							/*WARN USER to check data */
		 
				for (i=0;i<filelength;i++) {		
					
					if (     (*(input + ( i * COLS + col4))) == groupticker	  )		{			/*same group*/
						(*(input + ( i * COLS + col1)))  =   ( (  (*(input + ( i * COLS + col1))) + 2) * (-1));  						 
					}
				}
			}
 
		} /* loop for chossing best V */


		else {
			printf("Q");
		}
 
		if (groupnumber <= 0)	{
			logicerror(1);
		}

		groupticker = groupticker + 1;		

	}  /*loop for every group */




	/* Part Three -- OUTPUT the best in each group to file -- print output to file "predWORK" *****/

	sprintf(fullname, "%s.txt", outfile);
	/** OPEN predWORK FILE for writing **/

	output = fopen(fullname, "a+"); 
	if (output == NULL) 	 {
		printf("\nFailure to Create PredWORK File.\n");
		exit(2);
	}


	for (i=0;i<filelength;i++) {


		/**** for debugging*/
		printf("\n%d\t%d\t%d\t%.5f\t%.5f\t%.1f", 
				(int) (*(input + (i * COLS + 0))),
				(int) (*(input + (i * COLS + 1))), 
				(int) (*(input + (i * COLS + 2))), 
				(float) (*(input + (i * COLS + 5))),
				(float) (*(input + (i * COLS + 7))),
				(float) (*(input + (i * COLS + 6)))  ); 
		/***** end of debugging screen */


		if (  (*(input + (i * COLS + col1)))  != NEGONE)	{

				fprintf(output, "%d\t%d\t%.3f\t%.1f\t%.0f\t%.3f\n", /*see definitions.h*/ 
					(int) (*(input + (i * COLS + 0))),
					(int) (*(input + (i * COLS + 1))), 
					(float) (*(input + (i * COLS + col1))),			/* V */
					(float) (*(input + (i * COLS + col2))),			/* -logP */
					(float) (*(input + (i * COLS + col3))),			/* df */
					(float) (*(input + (i * COLS + col5)))    );	/* sens */
				
				printcounter = printcounter + 1;
		}
	}
  
	fclose(output);

}




/************************************************/
/************  PATHFINDER  **********************/
/************************************************/
 



void pathfinder(int filepos, int group, double *datafile, int filemax)	{

/*X = value of first position; Y= value of second position, filemax = total number of items in PVSDscan file*/


/*COLS is global offset given file shape  */
/* Need to put data somewhere */
/* Need to mark old data -- eliminator2 function */

	void pathfinder(int filepos, int group, double *datafile, int filemax); /*since recurisive, function prototype*/

	int xup;
	int xdown;
	int yup;
	int ydown;

	xup = filepos;
	xdown = filepos;
	yup = filepos;
	ydown = filepos;


	(*(datafile + (filepos * COLS + col4))) = group;
	


	/* Four loops for four 'directions */

	/* Loop One */

	while (xup < (filemax-1))	{
		
		/* if X finds a match seeking up in either position 0 or position 1*/

		if (   
			     (    (*(datafile + (filepos * COLS + 0))) == (*(datafile + ( (xup+1) * COLS + 0)))  
			       || (*(datafile + (filepos * COLS + 0))) == (*(datafile + ( (xup+1) * COLS + 1)))    
																								   )
			 &&  (	(*(datafile + ( (xup+1) * COLS + col4))) >= 0									   )
					
		   )	{


			pathfinder( (xup+1), group, datafile, filemax);

		}

		xup +=1;					 

	}



	/*Loop Two */


	while (xdown > 0)	{
		
		/* if X finds a match seeking down in either position 0 or position 1*/
 

		if (   
			     (    (*(datafile + (filepos * COLS + 0))) == (*(datafile + ( (xdown-1) * COLS + 0)))  
			       || (*(datafile + (filepos * COLS + 0))) == (*(datafile + ( (xdown-1) * COLS + 1)))    
																								   )
			 &&  (	(*(datafile + ( (xdown-1) * COLS + col4))) >= 0									   )
					
		   )	{



			pathfinder( (xdown-1), group, datafile, filemax);

		}

		xdown -=1;					 

	}


	/* Loop three */

	while (yup < (filemax-1))	{
		
		/* if y finds a match seeking up in either position 0 or position 1*/

		if (   
			     (    (*(datafile + (filepos * COLS + 1))) == (*(datafile + ( (yup+1) * COLS + 0)))  
			       || (*(datafile + (filepos * COLS + 1))) == (*(datafile + ( (yup+1) * COLS + 1)))    
																								   )
			 &&  (	(*(datafile + ( (yup+1) * COLS + col4))) >= 0									   )
					
		   )	{


			pathfinder( (yup+1), group, datafile, filemax);

		}

		yup +=1;					 

	}



	/*Loop four */


	while (ydown > 0)	{
		
		/* if y finds a match seeking down in either position 0 or position 1*/
 

		if (   
			     (    (*(datafile + (filepos * COLS + 1))) == (*(datafile + ( (ydown-1) * COLS + 0)))  
			       || (*(datafile + (filepos * COLS + 1))) == (*(datafile + ( (ydown-1) * COLS + 1)))    
																								   )
			 &&  (	(*(datafile + ( (ydown-1) * COLS + col4))) >= 0									   )
					
		   )	{



			pathfinder( (ydown-1), group, datafile, filemax);

		}

		ydown -=1;					 

	}
	 
}



/************************************************/
/************  LOGICERROR  **********************/
/************************************************/

void logicerror(int code)	{

	void	exitprogram(); 

	printf("\n\n\nPROGRAMMING LOGIC ERROR !!!!!! \n\n\n");
	printf("\nContact phillip.pang@stanfordalumni.org");
	printf("\n\nReference the following error:");

	if (code == 1)	{
		printf("\ngroupnumber <= 0; chain elimination flaw.");
	}

	exitprogram();
	exit(34);
}
