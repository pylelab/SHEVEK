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

 
int		**AR_analysis(char *input, int offset, int numrows, int column1, int column2, int *num_seq);
void	arcalc(int **nn, int ni, int nj, float *chisq, float *df, double *prob, float *cramrv, float *ccc, float **ar);
void	apply_thresholds(void);
int		cell_finder(char *input,int offset, int numrows, int num_start, char *unique1, char *unique2, int count1, int count2, 
				 int column1, int column2,int **chi_matx2,float **ar_matrix);
int		**chi_analysis(char *input, int start, int stop);
int		cochrantest(int **chi_matrix, int *rowtot, int *coltot, float *expctd, int numrows, int numcols);
float	**convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
void	crosstab(char *input, int **chi_matrix, char *aa1, char *aa2, int var1, int var2, int count1, int count2);
void	crosstab2(char *input, int **chi_matrix, char *aa1, char *aa2, int var1, int var2, int count1, int count2);
void	distribution_analyzer(float *VThresh, float *PThresh, float *DFThresh);
void	eliminator(double *input,int filelength);
void	eliminator2(double *input,int filelength);
double	estExact(int *ROWmatrix, double *dmatrix, double PRECISION, int numrows, int numcols, int *Echeck);
void	exitprogram();
int		fexactmain (int nrow, int ncol, double *table, double *Eprob);
double	fastexp2 (double value);
int		find_unique_elements(char *input, char *aa, int *gap, int var);
int		find_unique_elements2(char *input, int numrows, int offset, char *aa, int *gap, int var);
int		find_unique_seq(int *inp, int *out, int ns);
void	gap_rectifer(int **chi_matrix,int gap1, int gap2,int *cnt1, int *cnt2);
void 	gap_rectifer1(int **chi_matrix,int gap1, int *cnt1, int *cnt2, char *u1);
void 	gap_rectifer2(int **chi_matrix,int gap2, int *cnt1, int *cnt2, char *u2);
float	GetFloatTimer (float timer);
void	intro();
int	main ();
void	matrixconverter(int **chi_matrix, int *matrix, double *dmatrix, int count1, int count2, int flag);
void	misalign_identifier(char *input, int offset, int numrows);
void	misalignment_manager(char *input,int offset,int numrows);
void	message(int number);
char	*openfile();
double	*openfile2(int pass, int filelength);
double	*openfile3(int pass,float *actpmax);
double	*openfile4();
void	pathfinder(int fileposition, int group, double *datafile, int filemax);
void	positionrelater(char *input);
int		prelimscan(double *input, int filelength);
double	rcont2f (void);
int		rcont2s (int nrow, int ncol,int *_nrowt,int *_ncolt,int *_matrix);
void	read_inout(FILE *ifile);
double	*read_input3(FILE *ifile, int pass, float *actpmax);
double	*read_input4(FILE *ifile);
void	redunelim(void);
int		sequence_finder(char *input, int offset, int numrows, int start, char char1, char char2, int column1, int column2);
void	sorter(float *A, int left, int right);
void	score_manager(char *input);
void	screener(float *VThresh, float *PThresh, float *DFThresh);
void	swapper (float *A, int i, int j);

double	*dvector(long nl, long nh);
void	free_dvector(double *v, long nl, long nh);
void	exitprogram(); 
float   chipang(int **nn, int ni, int nj, float *chisq, float *df, double *prob, double *cramrv, float *ccc);
void	file_manager(void);
void	logicerror(int code);
char	*read_input(FILE *ifile);
static int pcomp(const void *p1, const void *p2);
static int vcomp(const void *p1, const void *p2);
static int dfcomp(const void *p1, const void *p2);
/*float	sensetest(int **nn, int ni, int nj); */
/*int	chisquare(int **nn, int ni, int nj, float *chisq, float *df, double *prob, double *cramrv, float *ccc); */




