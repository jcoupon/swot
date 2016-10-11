/*
 *    utils.h
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */

#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#include "mpi.h"
#include <gsl/gsl_integration.h>
#include "fitsio.h"

/* 	useful constants */
#define PI    3.14159265358979323846
#define TWOPI 6.283185307179586476925287
#define c     299792.458
#define G     6.67300e-11

/* 	WMAP5 Cosmology */
#define H0 72 /* all values definied with h = 0.72 */
#define Omega_M 0.258
#define Omega_L 0.742

/*		this set of utils is for file reading */
#define NFIELD 1000
#define NCHAR  100
#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)     atoi(array+NCHAR*(col-1))
#define getCharValue(array,col)    array+NCHAR*(col-1)
#define getLine(array,i)           array+NFIELD*NCHAR*i

/* 	global definitions */
#define FAILURE     0
#define SUCCESS     1
#define ODD         0
#define EVEN        1
#define COMO        0
#define THETA       1
#define PHYS        2

#define AUTO        0
#define CROSS       1
#define AUTO_WP     2
#define CROSS_WP    3
#define GGLENS      4
#define AUTO_3D     5
#define CROSS_3D    6
#define NUMBER      7

#define LS          0
#define NAT         1
#define HAM         2
#define PEEBLES     3
#define MASTER      0
#define FIRSTCALL   1
#define ROOT        0
#define JACKKNIFE   0
#define BOOTSTRAP   1
#define SUBSAMPLE   2

/* 	maximum stack "level" permitted. This means
 * 	that up to 2^NSTACK objects can be treated */
#define NSTACK 64

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN
#define SQUARE(a) ((a)*(a))


/* 	tree structure */
#define NLEAF        1 /* DO NOT CHANGE THIS */
#define NODE_STRIDE  2


/* 	BEWARE: avoid global variables as much as possible */
int NDIM;
char MYNAME[100];
double EPS;


/*
 *		new structures and global variables
 */


#define NIDSMAX   7 /* maximum number of columns */

/* 	data point structure */
typedef struct Point
{
  /* 	number of points */
  long N;

  /* 	dimension along which they are sorted
     	(ascending order). -1 if not.  */
  int dim;

  /* 	coordinates: Point.x[NDIM*n+DIM].
   * 	For gglensing and wp(rp), redshift is DIM = 3 */
  double *x;

  /* only used for gg lensing, Point.*[n] */
  double *zerr; /* redshift error         */
  double *e1;   /* ellipticity            */
  double *e2;
  double *w;    /* ellipticity weight     */
} Point;


/* 	tree structure */
typedef struct Tree
{
  /* 	total number of nodes */
  long size;

  /* 	each of these have "size" elements (2 times for cosx and sinx) */
  long *left, *right;
  double *cosx, *sinx, *r, *distComo;

  /* 	total number of objects per  sample
   * 	(for the whole sample it is tree.N[ROOT])
   */
  double *N, *Ntot;

  /* 	statistical weights for each sample */
  unsigned char *w;

  /* 	point with mean properties of the node (pointer to
   * 	a series of "size" elements arrays) */
  Point point;

} Tree;

typedef struct Mask
{
  /*		has NDIM*nsamples size */
  double *min, *max;

  /*		weights - has nsamples size */
  unsigned char  *w;
} Mask;

/* 	result structure (corr. function and errors) ----------- */
typedef struct Result
{
  /*		total number of points for the normalization,
   *		size = (nsamples+1)
   */
  double *N1, *N2;

  /* 		number of pairs, size = nbins*(nsamples+1) or nbins*nbins*(nsamples+1).
   * 		NN_s is for xi(s)
   */
  double *NN, *NN_s;

  /* 		for GG lensing, size = (nbins*nsamples+1) */
  double *GG, *w;

  /*		for GG lensing, size = nbins */
  double *Nsources, *e2, *meanR;

} Result;

/* 	config structure */
typedef struct Config
{
	/* 		correlation function options.
	 * 		Default proj(ection) for gg lensing is
	 * 		comoving (COMO), for auto.cross is theta (THETA).
	 * 		Change to PHYS for physical coordinates*/
	int cov_mat, estimator, nbins, nbins_pi,
	 corr,  log, Ninfo, proj, xi, weighted, calib,  resample2D;
	double deltaz, min, max, Delta, Delta_pi, OA, pi_max;

	/* 		error method JACKKNIFE or BOOTSTRAP */
	int err, nsamples, nsub;

	/*	 	fits files */
	int fits;

	/* 		Print tree */
	int printTree, printTreeAndExit;

	/* 		seed */
	size_t seed;

	/* 		write samples in files*/
	int printSamples;

	/* 		cosmology */
	double a[4];

	/* 		for mpi */
	int rank, size, verbose;

	/* 		input and output files information */
	char fileInName1[1000],  fileInName2[1000];
	char fileRanName1[1000], fileRanName2[1000];
	char fileOutName[1000],  RRInFileName[1000], RROutFileName[1000];

	/* 		column ids for input files */
	// int data1Id[NIDSMAX], data2Id[NIDSMAX];
	// int  ran1Id[NIDSMAX], ran2Id[NIDSMAX];
	char *data1Id[NIDSMAX];
	char *data2Id[NIDSMAX];
	char *ran1Id[NIDSMAX];
	char *ran2Id[NIDSMAX];

	/* 		angular separation function distAng[Spher,Cart] */
	double (*distAng)(const Tree *a, const long *i, const Tree *b, const long *j);

} Config;


/*
 * 	Utils - Numeric
 */

double distComo(double z,  const double a[4]);
double drdz(double x, void * params);
FILE *fopenAndCheck(const char *fileName, char *mode, int verbose);
int getStrings(char *line, char *strings, char *delimit, long *N);
double determineMachineEpsilon();
long randInt(long N);
void printCount(const long count, const long total, const long step, int verbose);
long lowestPowerOfTwo(long n, long *pow);
void comment(const Config para, char *commentString);
int checkFileExt(const char *s1, const char *s2);


/*
 *   Point routines
 */

void splitData(const Config para, const Point data, int dim, Point *dataLeft, Point *dataRight);
void quickSort(const Config para, Point data, int dim, long start, long end);
Point createPoint(const Config para, long N);
void freePoint(const Config para, Point point);
void copyPoint(const Config para, Point a, long i, Point b, long j);
void copyPointAddress(const Config para, Point *a, Point b, long shift);
void getMeanPoint(const Config para, Point a, long i, Point point);
void swapPoint(const Config para, Point point,long i, long j);
Point readCat(const Config para, char *fileInName, char *id[NIDSMAX], int weighted);

double distAngPointCart(const Config *para, const Point *a, const long *i, const Point *b, const long *j);
double distAngPointSpher(const Config *para, const Point *a, const long *i, const Point *b, const long *j);
double distAngSpher(const Tree *a, const long *i, const Tree *b, const long *j);
double distAngCart(const Tree *a, const long *i, const Tree *b, const long *j);
double distAngPoint(const Config *para, const Point *a, const long *i, const Point *b, const long *j);
double dist3D(const Tree *a, const long *i, const Tree *b, const long *j);


/*
 *    MPI routines
 */

void freeResult(const Config para, Result result);
void comData(const Config para, Point *data, long Ncpu, int dim, int firstCall);
void comResult(const Config para, Result result, long Ncpu, int split);







#endif
