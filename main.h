/*-------------------------------------------------------------------------*
 *main.h for swot (Super W Of Theta)                                       *
 *Jean Coupon, Alexie Leauthaud (2012)                                     *
 *-------------------------------------------------------------------------*/

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

/* useful constants */
#define PI    3.14159265358979323846
#define TWOPI 6.283185307179586476925287
#define c     299792.458
#define G     6.67300e-11

/* WMAP5 Cosmology */
#define H0 72 /* all values definied with h = 0.72 */
#define Omega_M 0.258
#define Omega_L 0.742

/* this set of utils is for file reading */
#define NFIELD 1000
#define NCHAR  100
#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getIntValue(array,col)     atoi(array+NCHAR*(col-1))
#define getCharValue(array,col)    array+NCHAR*(col-1)
#define getLine(array,i)           array+NFIELD*NCHAR*i

/* global definitions */
#define FAILURE     0
#define SUCCESS     1
#define ODD         0
#define EVEN        1
#define RADEC       0
#define CART        1
#define CART3D      2
#define COMO        0
#define THETA       1
#define PHYS        2

#define AUTO        0
#define CROSS       1
#define AUTO_WP     2
#define CROSS_WP    3
#define GGLENS      4

#define LS          0
#define NAT         1
#define HAM         2
#define MASTER      0
#define FIRSTCALL   1
#define ROOT        0
#define JACKKNIFE   0
#define BOOTSTRAP   1

/* maximum stack "level" permitted. This means 
 * that up to 2^NSTACK objects can be treated */
#define NSTACK 64

#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN
#define SQUARE(a) ((a)*(a))

/*----------------------------------------------------------------*
 * new structures and global variables                            *
 *----------------------------------------------------------------*/

/* data point structure -----------------------------------------*/
typedef struct Point
{
  /* number of points */
  long N;
  
  /* dimension along which they are sorted 
     (ascending order). -1 if not.  */
  int dim;

  /* coordinates: Point.x[NDIM*n+DIM]. 
   * For gglensing and wp(rp), redshift is DIM = 3 */
  double *x; 
  
  /* only used for gg lensing, Point.*[n] */
  double *zerr; /* redshift error         */
  double *e1;   /* ellipticity            */
  double *e2;
  double *w;    /* ellipticity weight     */
} Point;

/* tree structure */
#define NLEAF        1 /* DO NOT CHANGE THIS */
#define NODE_STRIDE  2

/* tree structure -------------------------------------------*/
typedef struct Tree
{
  /* total number of nodes */
  long size;
  
  /* each of these have "size" elements (2 times for cosx and sinx) */
  long *left, *right, *N;
  double *cosx, *sinx, *r, *distComo;
  
  /* total number of objects per  sample
   * (for the whole sample it is tree.N[ROOT])
   */
  long *Ntot;
  
  /* statistical weights for each sample */ 
  char *w;
  
  /* point with mean properties of the node (pointer to 
   * a series of "size" elements arrays) */ 
  Point point;

} Tree;

typedef struct Mask
{
  /* has NDIM*nsamples size */
  double *min, *max;
  
  /* weights - has nsamples size */
  char *w;
} Mask;

/* result structure (corr. function and errors) ----------- */
typedef struct Result
{
  /* total number of points for the normalization, 
   * size = (nsamples+1)
   */
  long *N1, *N2;
  
  /* number of pairs, size = nbins*(nsamples+1) */
  long *NN;
  
  /* for GG lensing, size = (nbins*nsamples+1) */
  double *GG, *w;
 
  /* for GG lensing, size = nbins */
  double *Nsources, *e2, *meanR;

} Result;

/* BEWARE: avoid global variables as much as possible */
int NDIM;
char MYNAME[100];
double EPS;


#define NIDSMAX   7 /* maximum number of columns */

/* config structure ---------------------------------------------*/
typedef struct Config
{
  /* correlation function options.
   * Default proj(ection) for gg lensing is
   * comoving (COMO), for auto.cross is theta (THETA).
   * Change to PHYS for physical coordinates*/
  int cov_mat, estimator, nbins, 
    corr, coordType, log, Ninfo, proj, xi;
  double deltaz, min, max, Delta, OA, pi_max; 
  
  /* error method JACKKNIFE or BOOTSTRAP */
  int err, nsamples;
  
  /* cosmology */
  double a[4];
  
  /* for mpi */
  int rank, size, verbose;
  
  /* input and output files information */
  char fileInName1[1000],  fileInName2[1000]; 
  char fileRanName1[1000], fileRanName2[1000];
  char fileOutName[1000];
  
  /* column ids for input files */
  int data1Id[NIDSMAX], data2Id[NIDSMAX];
  int  ran1Id[NIDSMAX], ran2Id[NIDSMAX];
  
  /* angular separation function distAng[Spher,Cart] */
  double (*distAng)(const Tree *a, const long *i, const Tree *b, const long *j);
  
} Config;

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

void autoCorr(Config para);
void crossCorr(Config para);
void ggCorr(Config para);

double wTheta(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int l);
double wp(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int j, int l);
Result Npairs(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall);
Result Npairs3D(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall);
Result gg(const Config *para,const Tree *lens, const long i, const Tree *source, const long j, int firstCall);
void corrLensSource(const Config *para, const Tree *lens, long i, const Tree *source, long j, double deltaTheta, Result result);
void freeResult(const Config para, Result result);

/*----------------------------------------------------------------*
 *Tree routines                                                   *
 *----------------------------------------------------------------*/

Tree buildTree(const Config *para, Point *data, Mask *mask, int dim, int firstCall);
void resample(const Config *para, const Point *data, int dim, Mask *mask, int firstCall);
void freeMask(const Config para, Mask mask);
void freeTree(Config para, Tree tree);
long splitTree(const Config *para, const Tree *tree1, const long root, const long Ncpu, int firstCall);
void printTree(const Config para, char *fileOutName, const Tree tree, long i, long NLeaf, int firstcall);
long countNodes(long N, long NLeaf);
double distAngPointCart(const Config *para, const Point *a, const long *i, const Point *b, const long *j);
double distAngPointSpher(const Config *para, const Point *a, const long *i, const Point *b, const long *j);
double distAngCart(const Tree *a, const long *i, const Tree *b, const long *j);
double distAngPoint(const Config *para, const Point *a, const long *i, const Point *b, const long *j);
double dist3D(const Tree *a, const long *i, const Tree *b, const long *j);

/*----------------------------------------------------------------*
 * MPI routines                                                   *
 *----------------------------------------------------------------*/

void comData(const Config para, Point *data, long Ncpu, int dim, int firstCall);
void comResult(const Config para, Result result, long Ncpu, int split);

/*----------------------------------------------------------------*
 *Configuration                                                   *
 *----------------------------------------------------------------*/

void initPara(int argc, char **argv, Config *para);
void setPara(char *field, char *arg, Config *para);
void checkArg(char *field, char *arg, Config *para);

/*----------------------------------------------------------------*
 *Utils - Numeric                                                 *
 *----------------------------------------------------------------*/

double distComo(double z,  const double a[4]);
double drdz(double x, void * params);
void splitData(const Config para, const Point data, int dim, Point *dataLeft, Point *dataRight);
void quickSort(const Config para, Point data, int dim, long start, long end);
Point createPoint(const Config para, long N);
void freePoint(const Config para, Point point);
void copyPoint(const Config para, Point a, long i, Point b, long j);
void copyPointAddress(const Config para, Point *a, Point b, long shift);
void getMeanPoint(const Config para, Point a, long i, Point point);
void swapPoint(const Config para, Point point,long i, long j);
Point readCat(const Config para, char *fileInName, int id[NIDSMAX]);
FILE *fopenAndCheck(const char *fileName, char *mode, int verbose);
int getStrings(char *line, char *strings, char *delimit, long *N);
double determineMachineEpsilon();
long randInt(long N);
void printCount(const long count, const long total, const long step, int verbose);
void comment(const Config para, char *commentString);
long lowestPowerOfTwo(long n, long *pow);

/*----------------------------------------------------------------*
 *obsolete                                                        *
 *----------------------------------------------------------------*/

Result Npairs2(const Config *para, const Tree *tree1, const long root1, const Tree *tree2, const long root2);
