#include <string.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "mpi.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_statistics.h>

#define FAILURE 0
#define SUCCESS 1


//Cosmology all value definied with h = 0.72
#define H0 72
#define Omega_M 0.258
#define Omega_L 0.742
#define c 299792.458
#define G 6.67300e-11

/* Alexie
Omega b  :    0.0438000
Omega m  :     0.258000
Omega l  :     0.742000
spectral index     0.963000
H0 :      72.0000
Sigma_8  :     0.796000
Critical density, z=0 :   1.43875e+11
*/

#define PI    3.14159265358979323846
#define TWOPI 6.283185307179586476925287

#define INF    1.0e30
#define ODD    0
#define EVEN   1
#define LEAF   0
#define NODE   1
#define RADEC  0
#define CART   1
#define CROSS  0
#define AUTO   1
#define GGLENS 2
#define LS     0
#define NAT    1
#define HAM    2
#define COMO   0
#define THETA  1
#define PHYS   2


#define MAX(x,y) ((x) > (y)) ? (x) : (y)
#define MIN(x,y) ((x) < (y)) ? (x) : (y)
#define ABS(a) ((a) < 0 ? -(a) : (a))
#define PARITY(a) (a)%2 ? ODD : EVEN

#define NFIELD 1000
#define NCHAR 100
#define getDoubleValue(array,col)  atof(array+NCHAR*(col-1))
#define getFloatValue(array,col)   atof(array+NCHAR*(col-1))
#define getIntValue(array,col)     atoi(array+NCHAR*(col-1))
#define getCharValue(array,col)    array+NCHAR*(col-1)
#define getLine(array,i)           array+NFIELD*NCHAR*i

#define NIDS 7

typedef struct Point
{
  unsigned char SplitDim, Nboots, *boots;
  double *x, *sinx, *cosx;
  double e1, e2, w, sigz;
  //long index;
} Point;

typedef struct Node
{
  Point *point;
  short Npoint;
  double radius,Deltaz;
  void *root, *Left, *Right;
  unsigned char type, Nboots;
  long id, *N;
} Node;


typedef struct LightNode
{
  void *root, *Left, *Right;
  unsigned char type;
} LightNode;

typedef struct Config
{
  int cov_mat, nboots, estimator, nbins,corr, coordType, log;
  
  /*For gg lensing, default: comoving (COMO). 
    For auto/cross, default: theta (THETA).
    change to PHYS for physical coordinate*/
  int proj;
  
  double (*distAng)(Point *,Point *);
  double min, max, Delta, OA;
  char fileInName1[1000], fileInName2[1000];
  char fileRanName1[1000], fileRanName2[1000];
  int data1Id[NIDS], data2Id[NIDS], ran1Id[NIDS], ran2Id[NIDS];
  char fileOutName[1000];
  char fileCovOutName[1000];
  double a[4],sigz;
} Config;

short NDIM = 2;
char MYNAME[100];
long IDERR;
double EPS = 1e-30;

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

int autoCorr(Config para, int rank, int size, int verbose);
int crossCorr(Config para, int rank, int size, int verbose);
Node *readAndCreateTree(Config para,char *fileName, int ids[NIDS], int rank, long *N);
void writeLandySzalay(Config para, long *D1D2, long *D1R1, long *D2R2, long *R1R2, long Ndata1, long Nrandom1, long Ndata2, long Nrandom2);
void writeNatural(Config para, long *D1D2, long *R1R2, long Ndata1, long Nrandom1, long Ndata2, long Nrandom2);
void writeHamilton(Config para, long *D1D2, long *D1R1, long *D2R2, long *R1R2, long Ndata1, long Nrandom1, long Ndata2, long Nrandom2);

/*----------------------------------------------------------------*
 *Tree routines                                                   *
 *----------------------------------------------------------------*/

Node  **getNodesForCpu(Config *para, Node *data, int count, int firstCall);
long *Npairs(Config *para, Node *node1, Node *node2, int rank, int size,  int firstCall, int verbose);
Node *createNode(Config para, Point *data, long N, int SplitDim, int firstCall);
void gg(Config *para, Node *source, Node *lens, double *SigR, double *w, int rank, int size,  int firstCall, int verbose);
void corrLensSource(Config *para, Point lens, Point source, double d, double *SigR, double *weight);

/*----------------------------------------------------------------*
 *Configuration                                                   *
 *----------------------------------------------------------------*/

int readPara(int argc, char **argv, Config *para);

/*----------------------------------------------------------------*
 *Utils - Numeric                                                 *
 *----------------------------------------------------------------*/

double dComo(double z, double a[4]);
double drdz(double x, void * params);
int getCatSize(FILE *fileIn);
void readCat(FILE *fileIn, int id[NIDS], double *data);
double distAngSpher(Point *a, Point *b);
double distAngSpherLog(Point *a, Point *b);
double distAngCart(Point *a, Point *b);
double distAngCartLog(Point *a, Point *b);
int comparePoints(const void *a,const void *b);
Point *dataToPoint(Config para, double *data, int N);
void cpyPoints(Point *a, Point *b);
void free_Point(Point *point, long N);
Point *createPoint(long N, short Nboots);
void free_Node(Node *node);
int compareDoubles(const void *a,const void *b);
double determineMachineEpsilon();
long determineIntError();
FILE *fopenAndCheck(const char *fileName,char *mode);
int getStrings(char *line, char *strings, char *delimit, long *N);
int roundToNi(double a);
void printCount(const long *count, const long *total, int verbose);
long randInt(long N);

