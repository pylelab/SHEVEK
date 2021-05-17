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


/**************************************************************/
/************   SCREENER function       ***********************/
/**************************************************************/
 

void screener(float *VThresh, float *PThresh, float *DFThresh)
{
	void distribution_analyzer(float *VThresh, float *PThresh, float *DFThresh);
	void	apply_thresholds(void);
	void	message(int number);
 
	message(1);		/*print to screen scoring message*/

	distribution_analyzer(VThresh,PThresh, DFThresh);   /*caclulate IQR and thus V&P Thresholds */



	printf("\n*************************************************");

	printf("\n\nSUGGESTED Lower -Log(P) threshold: %.1f", *PThresh);
	printf("\n\nSUGGESTED Lower V threshold: %.3f", *VThresh);
	printf("\n\nSUGGESTED protein DF threshold: %.1f", *DFThresh);


	printf("\n\n*************************************************");
	 

	message(2);		/*print to screen threshold warning message*/
	apply_thresholds();							 
	message(3);		/*print to screen prediction message*/

}

 