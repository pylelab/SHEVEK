/**************************************************/
/********from Numerical Recipies in C  ************/
/**************************************************/

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
 

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#define NR_END 1
#define FREE_ARG char*


void nrerror(const char* error_text)
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

 
double *dvector(long nl, long nh)
/* allocate a DOUBLE vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) {
		printf("FAILURE\n");
		exit(3);
	}
	return v-nl+NR_END;
}
 

void free_dvector(double *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	nh = nh;
	free((FREE_ARG) (v+nl-NR_END));
}


float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror((const char *)("allocation failure in vector()"));
	return v-nl+NR_END;
}


int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror((const char *)("allocation failure in vector()"));
	return v-nl+NR_END;
}




char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	char *v;

	v=(char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(char)));
	if (!v) nrerror((const char *)("allocation failure in cvector()"));
	return v-nl+NR_END;
}


float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror((const char *)("allocation failure 1 in matrix()"));
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror((const char *)("allocation failure 2 in matrix()"));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror((const char *)("allocation failure 1 in imatrix()"));
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror((const char *)("allocation failure 2 in imatrix()"));
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


int **subimatrix(int **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	int **m;

	oldch = oldch;

	/* allocate array of pointers to rows */
	m=(int **) malloc((size_t) ((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror((const char *)("allocation failure in submatrix()"));
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}



float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror((const char *)("allocation failure in convert_matrix()"));
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}



void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	nh = nh;
	free((FREE_ARG) (v+nl-NR_END));
}


void free_ivector(int *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	nh = nh;
	free((FREE_ARG) (v+nl-NR_END));
}



void free_cvector(char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	nh = nh;
	free((FREE_ARG) (v+nl-NR_END));
}


void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	nch = nch; nrh = nrh;
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	nch = nch; nrh = nrh;
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_subimatrix(int **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	nch = nch; ncl = ncl; nrh = nrh;
	free((FREE_ARG) (b+nrl-NR_END));
}


void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	nch = nch; ncl = ncl; nrh = nrh;
	free((FREE_ARG) (b+nrl-NR_END));
}

