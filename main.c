#include "main.h"

/*-------------------------------------------------------------------------*
 *swot (Super W Of Theta)  mpi version
 *Jean Coupon, Alexie Leauthaud (2012)
 *
 *Program to compute two-point correlation functions based on 
 *"divide and conquer" algorithms, mainly, but not limited to:
 *- data storage in a binary tree, 
 *- approximation at large scale,
 *- parellelization.
 *Supports auto and cross correlations, and galaxy-galaxy lensing.
 * 
 *TO DO:
 * - option to take into account the East-West orientation
 * - measure the tangential shear
 * - Do the randoms one by one or all together at once? -> Still need 
 *   to decide this question. NOT URGENT
 * 
 *Contributions:
 *- the algorithm to compute the number of pairs from a tree is 
 *  based on, but slightly different to, Martin Kilbinger's Ahtena code:
 *  http://www2.iap.fr/users/kilbinge/athena/
 *- the gal-gal lensing algorithm is based on
 *  Alexie Leauthaud's code 
 *  (see  Leauthaud et al. (2010),  2010ApJ...709...97L).
 *
 *Versions:
 *
 *v 0.12 May 1st [Jean]
 * - xi(r) 3D added (option -coord CART3D)
 * - default parameter file (swot -d) improved
 *
 *v 0.11 April 4th [Jean]
 * - data2 for sources catalogue
 * - more details in swot -d, -o now works
 * - no approx.: -OA no
 *
 *v 0.1 April 3rd [Jean]
 *- version from the old - memory monster - one.
 *  For previous versions, see main.c in legacy/
 *  option "nboots" has been replaced by "nsamples"
 *  and must be a power of two
 *-------------------------------------------------------------------------*/

int main(argc,argv)
     int argc;
     char **argv;
{

  int rank = 0, size = 1;
  /* MPI initialization:
   * "rank" is the id of the current cpu in use [0: master, size-1: last cpu]
   * "size" is the total number of cpus in use ("-np" option when calling mpirun) 
   */ 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  /* global configuration -> read config file and command line
   * note: all cpus read set para on their side because
   * sending structures is much more complex. 
   */
  Config para;
  para.size = size;
  para.rank = rank;
  initPara(argc,argv,&para);
  
  /* computation time begins */
  double t0 = MPI_Wtime();
  
  switch (para.corr){
  case AUTO:
    /* two-point autocorrelation function */
    autoCorr(para);
    break;
  case CROSS:
    /* two-point cross-correlation function */
    crossCorr(para);
    break;
  case GGLENS:
    /* galaxy-galaxy lensing two-point cross-correlation function */
    ggCorr(para);
    break;
  }
  
  /* computation time */
  MPI_Barrier(MPI_COMM_WORLD);
  if(para.verbose) fprintf(stderr, "%s: elapsed time: %5.2f s\n", MYNAME, MPI_Wtime() - t0);
  
  /* end of the main program */
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/

void autoCorr(Config para){
  /* Computes the auto correlation function.
   * for autocorrelation, each cpu correlates all nodes from the root 
   * with all nodes from a subnode (which is only a part of the tree).
   */
  
  int dimStart = 0;
  long i, j, k, l, n;
  Point data, random;
  Tree dataTree, randomTree;
  Result DD, DR, RR;
  Mask mask;
  
  /* read files */
  if(para.rank == MASTER){
    comment(para,"Reading fileRan1..."); random = readCat(para, para.fileRanName1, para.ran1Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random.N);}
    
    comment(para,"Reading fileIn1....");  data  = readCat(para, para.fileInName1, para.data1Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n",data.N);}
  }  
  
  /* resample, build masks */
  comment(para, "Resampling...");
  resample(&para, &random, dimStart, &mask, FIRSTCALL);
  
  /* send data */
  comment(para, "sending data...");  
  comData(para, &random, 0, dimStart, FIRSTCALL);
  comData(para, &data  , 0, dimStart, FIRSTCALL);
  
  /* grow trees */
  comment(para, "building trees...");
  randomTree = buildTree(&para, &random, &mask, dimStart, FIRSTCALL);   freePoint(para, random);
  dataTree   = buildTree(&para, &data, &mask, dimStart, FIRSTCALL);     freePoint(para, data);
  
  comment(para, "done.\n");
  
  /* divide and conquer */
  long nodeSlaveRan  = splitTree(&para, &randomTree, ROOT, para.size, FIRSTCALL);
  long nodeSlaveData = splitTree(&para, &dataTree,   ROOT, para.size, FIRSTCALL);

  /* DEBUGGING
  if(para.verbose){
    para.rank = 3;
    nodeSlaveData  = splitTree(&para, &dataTree, ROOT, para.size, FIRSTCALL);
    printTree(para, para.fileOutName, dataTree, nodeSlaveData, 1, FIRSTCALL);
    exit(-1);
  }*/
  
  /* compute pairs */
  comment(para, "RR...       "); RR = Npairs(&para, &randomTree, ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
  comment(para, "DR...       "); DR = Npairs(&para, &dataTree,   ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
  comment(para, "DD...       "); DD = Npairs(&para, &dataTree,   ROOT, &dataTree,   nodeSlaveData, FIRSTCALL);
    
  freeMask(para, mask);
  freeTree(para, randomTree);
  freeTree(para, dataTree);
  
  /* each slave sends the result and master sums everything up */
  comResult(para, RR, para.size, 0);
  comResult(para, DR, para.size, 0);
  comResult(para, DD, para.size, 0);
  
  /* print out results */
  if(para.rank == MASTER){
    
    /* mean w(theta) and errors */
    double *wmean = (double *)malloc(para.nbins*sizeof(double));
    double *err_r = (double *)malloc(para.nbins*sizeof(double));
    double *err_p = (double *)malloc(para.nbins*sizeof(double));
    
    double norm;
    switch(para.err){
    case JACKKNIFE: norm = (double)(para.nsamples - 1)/(double)(para.nsamples); break;
    case BOOTSTRAP: norm = 1.0/(double)(para.nsamples - 1); break;
    }										  
    
    /* w(theta) errors */
    for(i=0;i<para.nbins;i++){
      wmean[i] = err_r[i] = 0.0;
      if(para.nsamples > 1){
	for(l=0;l<para.nsamples;l++) wmean[i] += wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1)/(double)para.nsamples;
	for(l=0;l<para.nsamples;l++) err_r[i] += SQUARE(wmean[i]-wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1));
	/* resampling error */
	err_r[i] = sqrt(norm*err_r[i]);
      }
      /* poisson error ~1/N */
      err_p[i] = ABS(1.0+wTheta(para, para.estimator, DD, RR, DR, DR, i, 0))*(1.0/sqrt((double)DD.NN[para.nbins*0+i]) + 1.0/sqrt((double)RR.NN[para.nbins*0+i]));
    }
    
    FILE *fileOut = fopen(para.fileOutName,"w");
    switch(para.estimator){
    case LS:  fprintf(fileOut, "#Auto-correlation. Landy & Szalay estimator "); break;
    case NAT: fprintf(fileOut, "#Auto-correlation. Natural estimator ");        break;
    case HAM: fprintf(fileOut, "#Auto-correlation. Hamilton estimator ");       break;
    }
    switch(para.coordType){
    case RADEC: fprintf(fileOut, "(angular coordinates)\n"); break;
    case CART:  fprintf(fileOut, "(cartesian coordinates)"); break;
    }
    switch(para.err){
    case JACKKNIFE: fprintf(fileOut, "#Resampling: jackknife (%d samples)\n", para.nsamples); break;
    case BOOTSTRAP: fprintf(fileOut, "#Resampling: bootstrap (%d samples)\n", para.nsamples); break;
    }
    switch(para.proj){
    case THETA: fprintf(fileOut, "#  theta        w            err(total) err(resampling) err(poisson)"); break;
    case COMO:  fprintf(fileOut, "#  rp           w            err(total) err(resampling) err(poisson)"); break;
    }
    fprintf(fileOut, "       DD               DR            RR         Ndata      Nrandom\n");
    for(i=0;i<para.nbins;i++){
      if(DD.NN[i] > 0 && RR.NN[i] > 0){
	fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %15zd %15zd %15zd %10zd %10zd\n", 
		exp(para.min+para.Delta*(double)i+para.Delta/2.0), 
		wTheta(para, para.estimator, DD, RR, DR, DR, i, 0),
		sqrt(SQUARE(err_r[i])+SQUARE(err_p[i])), err_r[i], err_p[i],
		DD.NN[i], DR.NN[i], RR.NN[i],  DD.N1[0],  RR.N1[0]);
      }else{
	fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %15zd %15zd %15zd %10zd %10zd\n", 
		exp(para.min+para.Delta*(double)i+para.Delta/2.0), 
		0.0, 0.0, 0.0, 0.0,
		DD.NN[i], DR.NN[i], RR.NN[i],  DD.N1[0],  RR.N1[0]);
      }
    }
    fclose(fileOut);
    
    /* covariance matrix */
    if(para.cov_mat){ 
      double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
      FILE *fileCovOut = fopen(para.fileCovOutName,"w");
      for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
      for(i=0;i<para.nbins;i++){
	for(j=0;j<para.nbins;j++){
	  for(l=0;l<para.nsamples;l++){
	    cov[para.nbins*i+j] += norm*(wmean[i]-wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1))*(wmean[j]-wTheta(para, para.estimator, DD, RR, DR, DR, j, l+1));
	  }
	  /* add poisson error to diagonal */
	  if(i==j) cov[para.nbins*i+j] += SQUARE(err_p[i]);
	  /* write outfile */
	  if(DD.NN[i] > 0 && RR.NN[i] > 0 && DD.NN[j] > 0 && RR.NN[j] > 0){
	    fprintf(fileCovOut,"%g ", cov[para.nbins*i+j]);
	  }else{
	    fprintf(fileCovOut,"%f ", 0.0);
	  }
	}
	fprintf(fileCovOut,"\n");
      }
      fclose(fileCovOut);
      free(cov);
    }
    
    free(wmean);
    free(err_r);
    free(err_p);
  }
  
  freeResult(para, RR);
  freeResult(para, DR);
  freeResult(para, DD);
  
  return;
}

void crossCorr(Config para){
  /* Computes the cross-correlation function.
   * for autocorrelation, each cpu correlates all nodes from the root 
   * with all nodes from a subnode (which is only a part of the tree).
   */
  
  int swapped = 0, dimStart = 0;
  long i, j, k, l, n;
  Point data1, data2, random1, random2, tmp;
  Tree dataTree1, dataTree2, randomTree1, randomTree2;
  Result D1D2, D1R1, D2R2, R1R2;
  Mask mask1, mask2;
  
  /* read files */
  if(para.rank == MASTER){
    comment(para,"Reading fileRan1..."); random1 = readCat(para, para.fileRanName1, para.ran1Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random1.N);}
    
    comment(para,"Reading fileRan2..."); random2 = readCat(para, para.fileRanName2, para.ran2Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random2.N);}
    
    comment(para,"Reading fileIn1..."); data1 = readCat(para, para.fileInName1, para.data1Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n",data1.N);}
    
    comment(para,"Reading fileIn2..."); data2 = readCat(para, para.fileInName2, para.data2Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n",data2.N);}
    
    /* swap "1" and "2" if "1" is larger (to improve memory management)  */
    if(random1.N > random2.N){
      comment(para, "ATTENTION: fileRan1 is larger than fileRan2. Swapping \"1\" and \"2\" to save memory...\n");
      tmp = data1;   data1   = data2;   data2   = tmp;
      tmp = random1; random1 = random2; random2 = tmp;
      swapped = 1;
    }
  }  
  
  /* resample, build masks */
  comment(para, "Resampling...");
  resample(&para, &random1, dimStart, &mask1, FIRSTCALL);
  resample(&para, &random2, dimStart, &mask2, FIRSTCALL);
  
  /* send data */
  comment(para, "sending data...");  
  comData(para, &random1, 0, dimStart, FIRSTCALL);
  comData(para, &random2, para.size, dimStart, FIRSTCALL); /* 2 is partitioned among cpus */
  comData(para, &data1, 0, dimStart, FIRSTCALL);
  comData(para, &data2, 0, dimStart, FIRSTCALL);
  
  /* grow trees */
  comment(para, "building trees...");
  randomTree1 = buildTree(&para, &random1, &mask1, dimStart, FIRSTCALL); freePoint(para, random1);
  randomTree2 = buildTree(&para, &random2, &mask2, dimStart, FIRSTCALL); freePoint(para, random2);
  dataTree1 = buildTree(&para, &data1, &mask1, dimStart, FIRSTCALL);     freePoint(para, data1);
  dataTree2 = buildTree(&para, &data2, &mask2, dimStart, FIRSTCALL);     freePoint(para, data2);
  
  comment(para, "done.\n");                    
  
  /* divide and conquer */
  long nodeSlaveRan1  = splitTree(&para, &randomTree1, ROOT, para.size, FIRSTCALL);
  long nodeSlaveData2 = splitTree(&para, &dataTree2,   ROOT, para.size, FIRSTCALL);
  
  /* compute pairs */
  comment(para, "R1R2...       "); R1R2 = Npairs(&para, &randomTree1, ROOT, &randomTree2, ROOT,  FIRSTCALL);
  comment(para, "D1R1...       "); D1R1 = Npairs(&para, &dataTree1,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
  comment(para, "D2R2...       "); D2R2 = Npairs(&para, &dataTree2,   ROOT, &randomTree2, ROOT,  FIRSTCALL);
  comment(para, "D1D2...       "); D1D2 = Npairs(&para, &dataTree1,   ROOT, &dataTree2,   nodeSlaveData2, FIRSTCALL);
  
  freeMask(para, mask1);
  freeMask(para, mask2);
  freeTree(para, randomTree1);
  freeTree(para, randomTree2);
  freeTree(para, dataTree1);
  freeTree(para, dataTree2);
  
  /* each slave sends the result and master sums everything up */
  comResult(para, R1R2, para.size, 1); /* "1" to tell MASTER to sum up the total number         */
  comResult(para, D1R1, para.size, 0); /* of objects since R2 has been partitionned among cpus  */
  comResult(para, D2R2, para.size, 1);
  comResult(para, D1D2, para.size, 0);
  
  /* print out results */
  if(para.rank == MASTER){
    
    /* mean w(theta) and errors */
    double *wmean = (double *)malloc(para.nbins*sizeof(double));
    double *err_r = (double *)malloc(para.nbins*sizeof(double));
    double *err_p = (double *)malloc(para.nbins*sizeof(double));
    
    double norm;
    switch(para.err){
    case JACKKNIFE: norm = (double)(para.nsamples - 1)/(double)(para.nsamples); break;
    case BOOTSTRAP: norm = 1.0/(double)(para.nsamples - 1); break;
    }										  
    
    /* w(theta) errors */
    for(i=0;i<para.nbins;i++){
      wmean[i] = err_r[i] = 0.0;
      if(para.nsamples > 1){
	for(l=0;l<para.nsamples;l++) wmean[i] += wTheta(para, para.estimator, D1D2, R1R2, D1R1, D2R2, i, l+1)/(double)para.nsamples;
	for(l=0;l<para.nsamples;l++) err_r[i] += SQUARE(wmean[i]-wTheta(para, para.estimator, D1D2, R1R2, D1R1, D2R2, i, l+1));
	/* resampling error */
	err_r[i] = sqrt(norm*err_r[i]);
      }
      /* poisson error ~1/N */
      err_p[i] = ABS(1.0+wTheta(para, para.estimator, D1D2, R1R2, D1R1, D2R2, i, 0))*(1.0/sqrt((double)D1D2.NN[para.nbins*0+i]) + 1.0/sqrt((double)R1R2.NN[para.nbins*0+i]));
    }
    
    FILE *fileOut = fopen(para.fileOutName,"w");
    if(swapped) fprintf(fileOut, "#ATTENTION: \"1\" and \"2\" have been swapped to save memory.\n");
    switch(para.estimator){
    case LS:  fprintf(fileOut, "#Cross-correlation. Landy & Szalay estimator "); break;
    case NAT: fprintf(fileOut, "#Cross-correlation. Natural estimator ");        break;
    case HAM: fprintf(fileOut, "#Cross-correlation. Hamilton estimator ");       break;
    }
    switch(para.coordType){
    case RADEC: fprintf(fileOut, "(angular coordinates)\n"); break;
    case CART:  fprintf(fileOut, "(cartesian coordinates)"); break;
    }
    switch(para.err){
    case JACKKNIFE: fprintf(fileOut, "#Resampling: jackknife (%d samples)\n", para.nsamples); break;
    case BOOTSTRAP: fprintf(fileOut, "#Resampling: bootstrap (%d samples)\n", para.nsamples); break;
    }
    switch(para.proj){
    case THETA: fprintf(fileOut, "#  theta        w            err(total) err(resampling) err(poisson)"); break;
    case COMO:  fprintf(fileOut, "#  rp           w            err(total) err(resampling) err(poisson)"); break;
    }
    fprintf(fileOut, "       D1D2             D1R1             D2R2          R1R2       Ndata1     Nrandom1   Ndata2     Nrandom2\n");
    for(i=0;i<para.nbins;i++){
      if(D1D2.NN[i] > 0 && R1R2.NN[i] > 0){
	fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %15zd %15zd  %15zd %15zd %10zd %10zd %10zd %10zd\n", 
		exp(para.min+para.Delta*(double)i+para.Delta/2.0), 
		wTheta(para, para.estimator, D1D2, R1R2, D1R1, D2R2, i, 0),
		sqrt(SQUARE(err_r[i])+SQUARE(err_p[i])), err_r[i], err_p[i],
		D1D2.NN[i], D1R1.NN[i], D2R2.NN[i], R1R2.NN[i],  D1D2.N1[0],  R1R2.N1[0],  D1D2.N2[0],  R1R2.N2[0]);
      }else{
	fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %15zd %15zd  %15zd %15zd %10zd %10zd %10zd %10zd\n", 
		exp(para.min+para.Delta*(double)i+para.Delta/2.0), 
		0.0, 0.0, 0.0, 0.0,
		D1D2.NN[i], D1R1.NN[i], D2R2.NN[i], R1R2.NN[i],  D1D2.N1[0],  R1R2.N1[0],  D1D2.N2[0],  R1R2.N2[0]);
      }
    }
    fclose(fileOut);
    
    /* covariance matrix */
    if(para.cov_mat){ 
      double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
      FILE *fileCovOut = fopen(para.fileCovOutName,"w");
      for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
      for(i=0;i<para.nbins;i++){
	for(j=0;j<para.nbins;j++){
	  for(l=0;l<para.nsamples;l++){
	    cov[para.nbins*i+j] += norm*(wmean[i]-wTheta(para, para.estimator, D1D2, R1R2, D1R1, D2R2, i, l+1))*(wmean[j]-wTheta(para, para.estimator, D1D2, R1R2, D1R1, D2R2, j, l+1));
	  }
	  /* add poisson error to diagonal */
	  if(i==j) cov[para.nbins*i+j] += SQUARE(err_p[i]);
	  /* write outfile */
	  if(D1D2.NN[i] > 0 && R1R2.NN[i] > 0 && D1D2.NN[j] > 0 && R1R2.NN[j] > 0){
	    fprintf(fileCovOut,"%g ", cov[para.nbins*i+j]);
	  }else{
	    fprintf(fileCovOut,"%f ", 0.0);
	  }
	}
	fprintf(fileCovOut,"\n");
      }
      fclose(fileCovOut);
      free(cov);
    }
    
    free(wmean);
    free(err_r);
    free(err_p);
  }
  
  freeResult(para, R1R2);
  freeResult(para, D1R1);
  freeResult(para, D2R2);
  freeResult(para, D1D2);
  
  return;
}

void ggCorr(Config para){
  /* Computes the auto correlation function.
   * for autocorrelation, each cpu correlates all nodes from the root 
   * with all nodes from a subnode (which is only a part of the tree).
   */
  
  int dimStart = 0;
  long i, j, k, l, n;
  Point lens, source;
  Tree lensTree, sourceTree;
  Mask mask;
  Result result;
  
  /* read files */
  if(para.rank == MASTER){
    comment(para,"Reading sources file..."); source = readCat(para, para.fileInName2, para.data2Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd sources found).\n", source.N);}
    
    comment(para,"Reading lenses file...."); lens   = readCat(para, para.fileInName1, para.data1Id);
    if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd lenses found).\n", lens.N);}
  }  
  
  /* resampling = build masks TO DO: check that */
  comment(para, "Resampling...");
  resample(&para, &source, dimStart, &mask, FIRSTCALL);
  
  /* send data. Source catalogue is partitioned among cpus */
  comment(para, "sending data...");  

  comData(para, &source, para.size, dimStart, FIRSTCALL);
  comData(para, &lens  , 0, dimStart, FIRSTCALL);
  
  /* grow trees */
  comment(para, "building trees...");
  sourceTree = buildTree(&para, &source, &mask, dimStart, FIRSTCALL);   freePoint(para, source);
  lensTree   = buildTree(&para, &lens,   &mask, dimStart, FIRSTCALL);   freePoint(para, lens);
  
  comment(para, "done.\n");
  
  comment(para, "Correlating lenses with sources...       "); result = gg(&para, &lensTree, ROOT, &sourceTree, ROOT, FIRSTCALL);
  
  freeMask(para, mask);
  freeTree(para, sourceTree);
  freeTree(para, lensTree);
  
  /* each slave sends the result and master sums everything up */
  comResult(para, result, para.size, 0);
  
  /* print out results */
  if(para.rank == MASTER){
    
    /* resampling errors */
    double norm;
    switch(para.err){
    case JACKKNIFE: norm = (double)(para.nsamples - 1)/(double)(para.nsamples); break;
    case BOOTSTRAP: norm = 1.0/(double)(para.nsamples - 1); break;
    }										  
    
    double *GG_mean = (double *)malloc(para.nbins*sizeof(double));
    double *err_r   = (double *)malloc(para.nbins*sizeof(double));
    
    for(i=0;i<para.nbins;i++){
      err_r[i] = 0.0;
      GG_mean[i] = 0.0;
      if(para.nsamples > 1){
	for(l=0;l<para.nsamples;l++) GG_mean[i] += (-result.GG[para.nbins*(l+1)+i]/result.w[para.nbins*(l+1)+i])/(double)para.nsamples;
	for(l=0;l<para.nsamples;l++) err_r[i] += SQUARE(GG_mean[i] + result.GG[para.nbins*(l+1)+i]/result.w[para.nbins*(l+1)+i]);
	/* resampling error */
	err_r[i] = sqrt(norm*err_r[i]);
      }
    }
    
    double R, meanR;
    FILE *fileOut = fopen(para.fileOutName,"w");
    fprintf(fileOut, "#Gal-gal lensing. Sigma(R) vs R, linear approximation\n");
    switch(para.err){
    case JACKKNIFE: fprintf(fileOut, "#Resampling: jackknife (%d samples)\n", para.nsamples); break;
    case BOOTSTRAP: fprintf(fileOut, "#Resampling: bootstrap (%d samples)\n", para.nsamples); break;
    }
    switch(para.proj){
    case PHYS:  fprintf(fileOut, "#Coordinates: physical\n"); break;
    case COMO:  fprintf(fileOut, "#Coordinates: comoving\n"); break;
    }
    fprintf(fileOut, "#Cosmolgy: H0 = %g, Omega_M = %g, Omega_L = %g\n", para.a[0], para.a[1], para.a[2]);
    fprintf(fileOut, "#  R(Mpc)  SigR(Msun/pc^2) err(weights) err(resampling) Nsources       <R>          e2\n");
    for(i=0;i<para.nbins;i++){
      /* reminder: non-resampled value are stored from i=0 to nbins - 1 in result.w and result.GG */
      
      /* R and Rmean (weighted) */
      if(para.log){ 
	R     = exp(para.min+para.Delta*(double)i+para.Delta/2.0);
       	meanR = exp(result.meanR[i]/result.w[i]);
      }else{
	R     = para.min+para.Delta*(double)i+para.Delta/2.0;
	meanR = result.meanR[i]/result.w[i];
      }
      
      /* don't divide by zero if there are too few objects in the bin */
      if(result.Nsources[i] > 3.0){
	fprintf(fileOut,"%12.7f %12.7f %12.7f %12.7f %15zd %12.7f %12.7f\n", 
		R, -result.GG[i]/result.w[i],				\
		sqrt(1.0/result.w[i]), 	err_r[i],			\
		(long)result.Nsources[i],				\
		meanR,							\
		result.e2[i]/result.w[i]);
      }else{
	fprintf(fileOut,"%12.7f %12.7f %12.7f %12.7f %15zd %12.7f %12.7f\n", 
		R, 0.0,	0.0, 0.0, (long)result.Nsources[i], 0.0, 0.0);
      }
    }
    
    fclose(fileOut);
    
    /* covariance matrix */
    if(para.cov_mat){ 
      double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
      FILE *fileCovOut = fopen(para.fileCovOutName,"w");
      for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
      for(i=0;i<para.nbins;i++){
	for(j=0;j<para.nbins;j++){
	  for(l=0;l<para.nsamples;l++){
	    cov[para.nbins*i+j] += norm*(GG_mean[i] + result.GG[para.nbins*(l+1)+i]/result.w[para.nbins*(l+1)+i])*(GG_mean[j] + result.GG[para.nbins*(l+1)+j]/result.w[para.nbins*(l+1)+j]);
	  }
	  /* write outfile */
	  if(result.Nsources[i] > 3.0 && result.Nsources[j] > 3.0){
	    fprintf(fileCovOut,"%g ", cov[para.nbins*i+j]);
	  }else{
	    fprintf(fileCovOut,"%f ", 0.0);
	  }
	}
	fprintf(fileCovOut,"\n");
      }
      fclose(fileCovOut);
      free(cov);
    }
  }
  
  freeResult(para, result);
  
  return;
}


double wTheta(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int l){
  /* i is the bin index. l is sample 0 to 256 (0: no resampling, 1 -> 256: bootstrap or jackknife samples) ) */
  
  /* initialization */
  double Norm1 = (double)(R1R2.N1[l]*(R1R2.N2[l]-1))/(double)(D1D2.N1[l]*(D1D2.N2[l]-1));
  double Norm2 = (double)(R1R2.N2[l]-1)/(double)D1D2.N1[l];
  double Norm3 = (double)(R1R2.N1[l]-1)/(double)D1D2.N2[l];
  double Norm4 = (double)(D1D2.N2[l]*R1R2.N2[l])/(double)((R1R2.N2[l]-1)*(D1D2.N2[l]-1));
  
  double result;
  
  switch(estimator){
  case LS:  /* Landy and Szalay */
    result  =  Norm1*(double)D1D2.NN[para.nbins*l+i]/(double)R1R2.NN[para.nbins*l+i];
    result += -Norm2*(double)D1R1.NN[para.nbins*l+i]/(double)R1R2.NN[para.nbins*l+i];
    result += -Norm3*(double)D2R2.NN[para.nbins*l+i]/(double)R1R2.NN[para.nbins*l+i] + 1.0;
    break;
  case NAT: /* Natural */
    result = Norm1*(double)D1D2.NN[para.nbins*l+i]/(double)R1R2.NN[para.nbins*l+i] - 1.0;
    break;
  case HAM: /* Hamilton */
    result = Norm4*(double)D1D2.NN[para.nbins*l+i]*(double)R1R2.NN[para.nbins*l+i]/((double)D1R1.NN[para.nbins*l+i]*(double)D2R2.NN[para.nbins*l+i]) - 1.0;
    break;
  }
  
  return result;
}

#define leaf(tree,node) ((tree)->left[(node)] == -1)
#define node(tree,node) ((tree)->left[(node)] > -1)

Result Npairs(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall){ 
  /* returns the number of pairs if tree1 != tree2, and twice 
   * the number of pairs if tree1 == tree2, so that the treatment of pairs 
   * in w(theta) estimators is identical (numerical trick).
   */
  
  long NN, k = 0, l;
  static Result result;
  static long count, total;
  double deltaTheta;
  
  if(firstCall){
    count = 0;
    total = tree1->N[i]*tree2->N[j];
    result.NN = (long *)malloc(para->nbins*(para->nsamples+1)*sizeof(long));
    for(k = 0; k < para->nbins*(para->nsamples+1); k++) result.NN[k] = 0;
    
    /* number of points (used in wTheta(...) for the normalization) */
    result.N1 = (long *)malloc((para->nsamples+1)*sizeof(long));
    result.N2 = (long *)malloc((para->nsamples+1)*sizeof(long));
    result.N1[0] = tree1->N[ROOT];
    result.N2[0] = tree2->N[ROOT];
    for(l=0;l<para->nsamples;l++){
      result.N1[l+1] = tree1->Ntot[l];
      result.N2[l+1] = tree2->Ntot[l];
    }
  }
  
  if(tree1 == tree2 && i > j) return result;
  
  deltaTheta = para->distAng(tree1, &i, tree2, &j);
  
  if(node(tree1, i) && tree1->r[i]/deltaTheta > para->OA){
    if(node(tree2, j) && tree2->r[j]/deltaTheta > para->OA){
      Npairs(para, tree1, tree1->left[i],  tree2, tree2->left[j],  0); 
      Npairs(para, tree1, tree1->right[i], tree2, tree2->left[j],  0); 
      Npairs(para, tree1, tree1->left[i],  tree2, tree2->right[j], 0); 
      Npairs(para, tree1, tree1->right[i], tree2, tree2->right[j], 0); 
    }else{
      Npairs(para, tree1, tree1->left[i], tree2,  j,  0); 
      Npairs(para, tree1, tree1->right[i], tree2, j,  0);
    }
  }else if(node(tree2,j) && tree2->r[j]/deltaTheta > para->OA){
    Npairs(para, tree1, i, tree2, tree2->left[j],   0); 
    Npairs(para, tree1, i, tree2,  tree2->right[j], 0); 
  }else{
    if(para->log) deltaTheta = log(deltaTheta);
    k = floor((deltaTheta - para->min)/para->Delta);
    if(0 <= k && k < para->nbins){
      NN = tree1->N[i]*tree2->N[j];
      result.NN[k] += NN;
      for(l=0;l<para->nsamples;l++){
	result.NN[para->nbins*(l+1) + k] += NN*tree1->w[para->nsamples*i + l]*tree2->w[para->nsamples*j + l];
      }
    }
    count += tree1->N[i]*tree2->N[j];
    printCount(count,total,10000,para->verbose);
  }
  
  if(firstCall && para->verbose) fprintf(stderr, "\b\b\b\b\b\b\b%6.2f%%\n",100.0); 
  
  if(firstCall && tree1 == tree2){
    for(k = 0; k < para->nbins*(para->nsamples+1); k++) result.NN[k] *= 2;
  }
  
  return result;
}

Result gg(const Config *para, const Tree *lens, const long i, const Tree *source, const long j, int firstCall){
  /* Computes the galaxy-galaxy two-point correlation function. */
  
  long k;
  double deltaTheta;
  static Result result;
  static long total, count;
  
  if(firstCall){
    total = lens->N[i]*source->N[j];
    count = 0;
    result.GG = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));
    result.w  = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));
    for(k = 0; k < para->nbins*(para->nsamples+1); k++){
      result.GG[k] = 0.0;
      result.w[k]  = 0.0;
    }
    result.Nsources = (double *)malloc(para->nbins*sizeof(double));
    result.meanR    = (double *)malloc(para->nbins*sizeof(double));
    result.e2       = (double *)malloc(para->nbins*sizeof(double));
    for(k = 0; k < para->nbins; k++){
      result.Nsources[k] = 0.0; 
      result.meanR[k]    = 0.0;
      result.e2[k]       = 0.0;
    }
  }
  
  deltaTheta = para->distAng(lens, &i, source, &j);
  
  /* correlate if (leaf OR size/d < OA), otherwise go down the tree. Redshift is point.x[NDIM*i + 2] */
  if(node(lens, i) && lens->r[i]/deltaTheta > para->OA ) {
    if(node(source, j) && source->r[j]/deltaTheta > para->OA){
      gg(para, lens, lens->left[i],  source, source->left[j],  0);
      gg(para, lens, lens->right[i], source, source->left[j],  0);
      gg(para, lens, lens->left[i],  source, source->right[j], 0);
      gg(para, lens, lens->right[i], source, source->right[j], 0);
    }else{
      gg(para, lens, lens->left[i],  source, j, 0);
      gg(para, lens, lens->right[i], source, j, 0);
    }
  }else if(node(source, j) && source->r[j]/deltaTheta > para->OA){
    gg(para, lens, i, source, source->left[j],  0);
    gg(para, lens, i, source, source->right[j], 0);
  }else{
    corrLensSource(para, lens, i, source, j, deltaTheta, result);
    count += lens->N[i]*source->N[j];
    printCount(count, total, 10000, para->verbose);
  }
  
  if(firstCall && para->verbose) fprintf(stderr, "\b\b\b\b\b\b\b%6.2f%%\n",100.0); 
  
  return result;
}

#define RA_lens        ((lens)->point.x[NDIM*(i)+0])
#define DEC_lens       ((lens)->point.x[NDIM*(i)+1])
#define z_lens         ((lens)->point.x[NDIM*(i)+2])
#define zerr_lens      ((lens)->point.zerr[(i)])
#define distComo_lens  ((lens)->distComo[(i)])

#define RA_source        ((source)->point.x[NDIM*(j)+0])
#define DEC_source       ((source)->point.x[NDIM*(j)+1])
#define z_source         ((source)->point.x[NDIM*(j)+2])
#define zerr_source      ((source)->point.zerr[(j)])
#define e1_source        ((source)->point.e1[(j)])
#define e2_source        ((source)->point.e2[(j)])
#define w_source         ((source)->point.w[(j)])
#define distComo_source  ((source)->distComo[(j)])

void corrLensSource(const Config *para, const Tree *lens, long i,const Tree *source, long j, double deltaTheta, Result result){
  /* See Leauthaud et al. (2010),  2010ApJ...709...97L */
  
  double dA, dR, invScaleFac;
  long k, zero = 0;
  
  /* very quick tests? But might be time consuming, too. The key 
   * is to find the best balance...
   * e.g. if(z_source < z_lens) return;
   */
  
  switch(para->proj){
  case COMO:
    invScaleFac = 1.0 + z_lens;
    break;
  case PHYS:
    invScaleFac = 1.0;
    break;
  case THETA:
    printf("Wrong coordinate system for gal-gal lensing\n");
    exit(EXIT_FAILURE);
    break;
  }
  
  /* dA = phy_dis/angular_size_in_radians = tran_como_dis/(1+z)*/
  dA  = distComo_lens/(1.0 + z_lens);               /* Angular diameter distance in physical coordinates */
  dR  = deltaTheta*dA*PI/180.0;                     /* Transverse distance in phys coordinates (Mpc)     */
  dR *= invScaleFac;                                /* If coordinate in comoving system, divide by a     */
  
  if(para->log) dR = log(dR);
  k = floor((dR - para->min)/para->Delta);
  if(0 <= k && k < para->nbins && z_source > z_lens + zerr_lens + zerr_source && z_source > zerr_lens + para->deltaz){
  //DEBUGGING if(0 <= k && k < para->nbins && z_source > z_lens + 0.1){
    /* Point A --------------------------------- */
    Point A = createPoint(*para, 1);
    A.x[0]    = RA_source;
    A.x[1]    = DEC_lens;
    double AS = distAngPointSpher(para, &A, &zero, &source->point, &j);
    double AL = distAngPointSpher(para, &A, &zero, &lens->point,   &i);
    /* to get correct sign for phi_gg ---------- */
    if(RA_source  > RA_lens)  AL = -AL;
    if(DEC_source < DEC_lens) AS = -AS;
    double phi_gg     = atan2(AS,AL);
    double cos2phi_gg = cos(2.0*phi_gg);
    double sin2phi_gg = sin(2.0*phi_gg);
    double e1         =  e1_source*cos2phi_gg + e2_source*sin2phi_gg;
    double e2         = -e1_source*sin2phi_gg + e2_source*cos2phi_gg;
    double DOS        = distComo_source;                 /* Note this is a comoving los distance,         */
    double DOL        = distComo_lens;                   /* but redshift factors cancel out               */
    double DLS        = distComo_source - distComo_lens; /* Approx. Omega_k = 0                           */
    double SigCritInv = DOL*DLS/DOS/1000.0/invScaleFac;  /* 1/SigCrit in Gpc                              */
    SigCritInv       /= 1.663e3;                         /* see Narayan & Bartelman pge 10                */
    SigCritInv       *= invScaleFac*invScaleFac;         /* If coordinates in comoving system             */
    
    /* ATTENTION: source->w and lens->w are resampling weights (boostrap or jackknife) whereas
     * w below (and further result.w) is the lensing weight (i.e. shape measurement error 
     * if inverse variance estimate) 
     */
    double w   = SigCritInv*SigCritInv*w_source;
    double GG  = e1*w/SigCritInv;
    
    long l;
    result.GG[k] += GG;
    result.w[k]  += w;
    for(l=0;l<para->nsamples;l++){
      result.GG[para->nbins*(l+1) + k] += GG*lens->w[para->nsamples*i + l]*source->w[para->nsamples*j + l];
      result.w[para->nbins*(l+1) + k]  += w*lens->w[para->nsamples*i + l]*source->w[para->nsamples*j + l];
    }    
    /* keep track of info per bin */
    result.Nsources[k] += 1.0; 
    result.meanR[k]    += dR*w;
    result.e2[k]       += e2*w/SigCritInv;
    freePoint(*para, A);
  }
}

#undef RA_lens
#undef DEC_lens
#undef z_lens
#undef zerr_lens
#undef distComo_lens

#undef RA_source
#undef DEC_source
#undef z_source
#undef zerr_source
#undef e1_source
#undef e2_source
#undef w_source
#undef distComo_source

#undef leaf(tree,node) 
#undef node(tree,node)

  
void freeResult(const Config para, Result result){
  
  switch(para.corr){
  case AUTO: case CROSS:
    free(result.N1);
    free(result.N2);
    free(result.NN);
    break;
  case GGLENS:
    free(result.GG);
    free(result.w);
    free(result.Nsources);
    free(result.meanR);
    free(result.e2);
    break;
  }
  
  return;
}

/*----------------------------------------------------------------*
 *Tree routines                                                   *
 *----------------------------------------------------------------*/

Tree buildTree(const Config *para, Point *data, Mask *mask, int dim, int firstCall){
  
  static long index, count;
  static Tree result;
  
  double d, dMax;
  long i, j, n, l, local_index, maxPoint;
  
  if(firstCall){
    
    /* number of nodes */
    result.size = countNodes(data->N, NLEAF);
    
    /* allocate memory slots (see structure description in main.h)*/
    result.left   = (long *)malloc(result.size*sizeof(long)); 
    result.right  = (long *)malloc(result.size*sizeof(long)); 
    result.N      = (long *)malloc(result.size*sizeof(long)); 
    result.r      = (double *)malloc(result.size*sizeof(double));
    result.w      = (char *)malloc(result.size*para->nsamples*sizeof(char));
    result.Ntot   = (long *)malloc(para->nsamples*sizeof(long));
    result.point  = createPoint(*para, result.size);
    if(para->coordType == RADEC){
      result.cosx   =  (double *)malloc(2*result.size*sizeof(double)); 
      result.sinx   =  (double *)malloc(2*result.size*sizeof(double)); 
    }
    if(para->corr == GGLENS){
      result.distComo = (double *)malloc(result.size*sizeof(double));
    }
    
    count = index = ROOT;                                       /* index of the root (= 0) */
    for(i=0;i<result.size*para->nsamples;i++) result.w[i] = 0;  /* set weights to 0        */
    for(i=0;i<para->nsamples;i++) result.Ntot[i] = 0;           /* initialize Ntot         */
  }
  
  local_index = count++;           /* index of THIS node             */
  result.N[local_index] = data->N; /* number of points for THIS node */
  
  /* get mean properties of points inside the node <=> "weighted" center */
  getMeanPoint(*para, result.point, local_index, *data);
  
  /* set weights from mask */
  for(j=0;j<para->nsamples;j++){  /* loop over subsamples */
    n = 0;
    for(i=0;i<NDIM;i++){
      n += (mask->min[NDIM*j + i] < result.point.x[NDIM*local_index + i] && 
	    result.point.x[NDIM*local_index + i] < mask->max[NDIM*j + i] );
    }
    if(n == NDIM){  /* "resul.point" is in the subsample "j" */
      for(i=0;i<para->nsamples;i++){  /* loop over resamplings  */
	result.w[para->nsamples*local_index + i] = mask->w[para->nsamples*i+j];
      }
      break;
    }
  }
  
  /* pre-computed quantities */
  if(para->coordType == RADEC){
    result.cosx[2*local_index+0] = cos(result.point.x[NDIM*local_index+0]*PI/180.0);
    result.cosx[2*local_index+1] = cos(result.point.x[NDIM*local_index+1]*PI/180.0);   
    result.sinx[2*local_index+0] = sin(result.point.x[NDIM*local_index+0]*PI/180.0);
    result.sinx[2*local_index+1] = sin(result.point.x[NDIM*local_index+1]*PI/180.0);   
  }
  if(para->corr == GGLENS){
    result.distComo[local_index] = distComo(result.point.x[NDIM*local_index+2], para->a);
  }
  
  /* compute node radius "r" for the open-angle approximation. r is defined as
   * the distance between the weighted center and the most distant point (angular 
   * separation). For speed purpose, the test uses the cartesian approximation, but 
   * the distance is then estimated accuratly according to the coordinate system.
   *
   * For the sake of "code readability" we make the assumption that box size
   * estimated on the projected coordinates should be a fair approximation 
   * of the 3-D size.
   */
  dMax = 0.0; maxPoint = 0;
  for(n=0;n<data->N;n++){
    d = distAngPointCart(para, &(result.point), &local_index, data, &n);
    if(d > dMax){
      maxPoint = n;
      dMax     = d;
    }
  }
  if(para->coordType == RADEC){
    result.r[local_index] = distAngPointSpher(para, &(result.point), &local_index, data, &maxPoint);
  }else{
    result.r[local_index] = d;
  }
  
  if(data->N > NLEAF){   /* node ----------------------------- */
    
    Point dataLeft, dataRight;
    
    /* split the node into two along the "dim" coordinate */
    splitData(*para, *data, dim, &dataLeft, &dataRight);
    
    /* next splitting coordinate */
    dim++; if(dim > NDIM-1) dim = 0;
    
    buildTree(para, &dataLeft, mask, dim, 0);    /* left child */
    result.left[local_index] = index;
    
    buildTree(para, &dataRight, mask, dim, 0);   /* right child */
    result.right[local_index] = index;
    
  }else{                 /* leaf ------------------------------ */
    result.left[local_index]  = -1;
    result.right[local_index] = -1;
    
    /* keep track of the total number of objects per sample */
    for(i=0;i<para->nsamples;i++){  
      result.Ntot[i] += NLEAF*result.w[para->nsamples*local_index + i];
    }
  }
  
  /* return the index of this node (index is static) */
  index = local_index;
  
  return result;
}

void resample(const Config *para, const Point *data, int dim, Mask *mask, int firstCall){
  /* splits data into a number of sub samples, resample it and builds up a 
   * mask with weights. (TO DO) Then send the mask (if MASTER) 
   * or receive it (if SLAVE). The mask is used in buildTree() 
   * to attribute a weight to each node. */

  static long depth, depthSample, count;
  long i, j, l, rank;

  if(firstCall){
    
    /* check up sample values value and evaluate depthSample */
    if(para->nsamples == 0){
      return;
    }else if(para->nsamples > 256){
      if(para->verbose) 
	fprintf(stderr,"\n%s: **ERROR** nsamples must be <= 256. Exiting...\n", MYNAME);
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }else if(lowestPowerOfTwo(para->nsamples, &depthSample) != para->nsamples){
      if(para->verbose)
	fprintf(stderr,"\n%s: **ERROR** nsamples must be a power of 2. Exiting...\n", MYNAME);
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    
    /* initialisation */
    mask->min = (double *)malloc(NDIM*para->nsamples*sizeof(double));  
    mask->max = (double *)malloc(NDIM*para->nsamples*sizeof(double));  
    mask->w   = (char *)malloc(para->nsamples*para->nsamples*sizeof(char));  
    
    if(para->rank == MASTER){
      
      if(para->nsamples > data->N){
	if(para->verbose) 
	  fprintf(stderr,"\n%s: **ERROR** nsamples must be < Npoints. Exiting...\n", MYNAME);
	exit(EXIT_FAILURE);
      }
      
      for(i=0;i<para->nsamples*para->nsamples;i++) mask->w[i] = 0;
      for(i=0;i<para->nsamples;i++){
	for(j=0;j<para->nsamples;j++){
	  /* i: resample draw index
	   * j: subsample index
	   */
	  switch(para->err){
	  case JACKKNIFE: mask->w[para->nsamples*i + j] = !(i == j);                break;
	  case BOOTSTRAP: mask->w[para->nsamples*i + randInt(para->nsamples)] += 1; break;
	  }
	}
      }
      depth = 0; /* depth of the node. Number of parallel nodes is 2^depth */
      count = 0; /* counter for mask array                                 */
    }
  }
  
  if(para->rank == MASTER){
    /* one level down */
    depth++;
    
    if(depth <= depthSample){
      
      Point dataLeft, dataRight;
      
      /* split the node into two along the "dim" coordinate */
      splitData(*para, *data, dim, &dataLeft, &dataRight);
      
      /* next splitting coordinate */
      dim++; if(dim > NDIM-1) dim = 0;
      
      resample(para, &dataLeft,  dim, mask, 0);
      resample(para, &dataRight, dim, mask, 0);
      
    }else{
      
      /* compute limits of the mask */
      for(dim=0;dim<NDIM;dim++){
    	mask->min[NDIM*count + dim] = data->x[NDIM*0+dim];
    	mask->max[NDIM*count + dim] = data->x[NDIM*0+dim];
      }
      for(i=1;i<data->N;i++){
    	for(dim=0;dim<NDIM;dim++){
    	  mask->min[NDIM*count + dim] = MIN(mask->min[NDIM*count + dim], data->x[NDIM*i+dim]);
    	  mask->max[NDIM*count + dim] = MAX(mask->max[NDIM*count + dim], data->x[NDIM*i+dim]);
    	}
      }
      
      count++;
    }
    
    /* one level up */
    depth--;
  }
  
  if(firstCall){
    MPI_Bcast(mask->min, NDIM*para->nsamples, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(mask->max, NDIM*para->nsamples, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(mask->w,   para->nsamples*para->nsamples, MPI_CHAR, MASTER, MPI_COMM_WORLD);
  }
  
  return;
}

void freeMask(const Config para, Mask mask){
  
  if(para.nsamples > 0){
    free(mask.min);
    free(mask.max);
    free(mask.w);
  }
  
}

void printTree(const Config para, char *fileOutName, const Tree tree, long i, long NLeaf, int firstCall){
  static FILE *fileOut;
  int l, dim;
  
  if(firstCall) fileOut = fopenAndCheck(fileOutName,"w", para.verbose);
  //if(firstCall) fileOut = stdout;
  
  if(tree.N[i] > NLeaf){
    printTree(para, fileOutName, tree, tree.left[i], NLeaf, 0);
    printTree(para, fileOutName, tree, tree.right[i], NLeaf, 0);
  }else{
    for(dim=0;dim<NDIM;dim++)   fprintf(fileOut,"%f ", tree.point.x[NDIM*i+dim]);
    for(l=0;l<para.nsamples;l++)  fprintf(fileOut,"%d ", tree.w[para.nsamples*i + l] );
    fprintf(fileOut,"%d \n", para.rank);
  }
  
  if(firstCall)  MPI_Barrier(MPI_COMM_WORLD);
 
}

void freeTree(const Config para, Tree tree){
  /* Frees the tree "tree". */
  
  free(tree.left);
  free(tree.right);
  free(tree.N);
  free(tree.r);
  free(tree.w);
  freePoint(para, tree.point);

  return;
}

long splitTree(const Config *para, const Tree *tree1, const long root, const long Ncpu, int firstCall){ 
  /* returns the id of the node to be computed by rank n. */
  
  static long count, rootRank; 
  
  if(firstCall) count = rootRank = 0;
  
  if(Ncpu > 1){
    switch(PARITY(Ncpu)){
    case EVEN:
      splitTree(para, tree1, tree1->left[root],  Ncpu/2, 0);
      splitTree(para, tree1, tree1->right[root], Ncpu/2, 0);
      break;
    case ODD:
      splitTree(para, tree1, tree1->left[root],  (Ncpu+1)/2, 0);
      splitTree(para, tree1, tree1->right[root], (Ncpu-1)/2, 0);
      break;
    }
  }else{
    if(count++ == para->rank) rootRank = root ; 
  }
  
  return rootRank;
}

long countNodes(long N, long NLeaf){
  /* Simply counts the number of nodes in a tree
   * from a data set of N points. NLeaf is the 
   * maximum number of objects in a leaf.
   */
  
  long result = 0;

  if(N > NLeaf){
    switch(PARITY(N)){ 
    case EVEN:
      result += 2*countNodes(N/2, NLeaf);
      break;
    case ODD:
      result += countNodes((N+1)/2, NLeaf);
      result += countNodes((N-1)/2, NLeaf);
      break;
    }
  }else{
    return 1;
  }
  
  return result + 1;
}

double distAngPointCart(const Config *para, const Point *a, const long *i, const Point *b, const long *j){
  /* returns the angular distance between points 
   * a[i] and b[i]. Cartesian coordinates.
   */
  
  double d0 = (b->x[NDIM*(*j)+0] - a->x[NDIM*(*i)+0]);
  double d1 = (b->x[NDIM*(*j)+1] - a->x[NDIM*(*i)+1]);
  
  return sqrt(d0*d0 + d1*d1);
}

double distAngPointSpher(const Config *para, const Point *a, const long *i, const Point *b, const long *j){
  /* returns the angular distance between points 
   * a[i] and b[i]. Spherical coordinates.
   */
  
  double sin2_ra  = 0.5*(1.0 - cos(a->x[NDIM*(*i)+0]*PI/180.0)*cos(b->x[NDIM*(*j)+0]*PI/180.0)
			 - sin(a->x[NDIM*(*i)+0]*PI/180.0)*sin(b->x[NDIM*(*j)+0]*PI/180.0));
  double sin2_dec = 0.5*(1.0 - cos(a->x[NDIM*(*i)+1]*PI/180.0)*cos(b->x[NDIM*(*j)+1]*PI/180.0) 
			 - sin(a->x[NDIM*(*i)+1]*PI/180.0)*sin(b->x[NDIM*(*j)+1]*PI/180.0));
  
  return 2.0*asin(sqrt(MAX(EPS/100.0, sin2_dec + cos(a->x[NDIM*(*i)+1]*PI/180.0)
			   *cos(b->x[NDIM*(*j)+1]*PI/180.0)*sin2_ra)))*180.0/PI;
}

double distAngSpher(const Tree *a, const long *i, const Tree *b, const long *j){
  /*Returns the angular distance between nodes 
   * a[i] and b[i]. Spherical coordinates.
   */

  double sin2_ra  = 0.5*(1.0 - a->cosx[2*(*i)+0]*b->cosx[2*(*j)+0] - a->sinx[2*(*i)+0]*b->sinx[2*(*j)+0]);
  double sin2_dec = 0.5*(1.0 - a->cosx[2*(*i)+1]*b->cosx[2*(*j)+1] - a->sinx[2*(*i)+1]*b->sinx[2*(*j)+1]);
  
  return 2.0*asin(sqrt(MAX(EPS/100.0,sin2_dec + a->cosx[2*(*i)+1]*b->cosx[2*(*j)+1]*sin2_ra)))*180.0/PI;
}

double distAngCart(const Tree *a, const long *i, const Tree *b, const long *j){
/*Returns the angular distance between nodes
 * a[i] and b[j]. Cartesian coordinates.
 */
  
  double d0 = (b->point.x[NDIM*(*j)+0] - a->point.x[NDIM*(*i)+0]);
  double d1 = (b->point.x[NDIM*(*j)+1] - a->point.x[NDIM*(*i)+1]);
  
  return sqrt(d0*d0 + d1*d1);
}


double dist3D(const Tree *a, const long *i, const Tree *b, const long *j){
/*Returns the angular distance between nodes
 * a[i] and b[j]. Cartesian coordinates in 3D.
 */
  
  double d0 = (b->point.x[NDIM*(*j)+0] - a->point.x[NDIM*(*i)+0]);
  double d1 = (b->point.x[NDIM*(*j)+1] - a->point.x[NDIM*(*i)+1]);
  double d2 = (b->point.x[NDIM*(*j)+2] - a->point.x[NDIM*(*i)+2]);
  
  return sqrt(d0*d0 + d1*d1 + d2*d2);
}


/*----------------------------------------------------------------*
 * MPI routines                                                   *
 *----------------------------------------------------------------*/

#define BASE 200
void comData(const Config para, Point *data, long Ncpu, int dim, int firstCall){
  /* Partitions and sends the data recursively if para.rank == master. 
   * Otherwise receives the data. Then returns the final data */
  
  MPI_Status status;
  
  if(para.rank == MASTER){
  
    static int count, rank;
    static Point *dataParent;
    
    Point dataLeft, dataRight;
    long NLeft, NRight;
  
    if(firstCall){
      count      = 0;
      rank       = para.size - 1;
      dataParent = data;
    }
    count++;
    
    if(Ncpu > 1){
      splitData(para, *data, dim, &dataLeft, &dataRight);
      
      /* next splitting coordinate */
      dim++; if(dim > NDIM-1) dim = 0;
      
      switch(PARITY(Ncpu)){
      case EVEN:
	comData(para, &dataLeft,  Ncpu/2, dim, 0);
	comData(para, &dataRight, Ncpu/2, dim, 0);
	break;
      case ODD:
	comData(para, &dataLeft,  (Ncpu+1)/2, dim, 0);
	comData(para, &dataRight, (Ncpu-1)/2, dim, 0);
	break;
      }
    }else{
      if(rank != MASTER){
	MPI_Send(&data->N, 1, MPI_LONG, rank, BASE+0, MPI_COMM_WORLD);
	MPI_Send(&dim, 1, MPI_INT, rank, BASE+1, MPI_COMM_WORLD);
	MPI_Send(data->x,data->N*NDIM, MPI_DOUBLE, rank, BASE+2, MPI_COMM_WORLD);
	
	if(para.corr == GGLENS){
	  MPI_Send(data->zerr, data->N, MPI_DOUBLE, rank, BASE+3, MPI_COMM_WORLD);
	  MPI_Send(data->e1  , data->N, MPI_DOUBLE, rank, BASE+4, MPI_COMM_WORLD);
	  MPI_Send(data->e2  , data->N, MPI_DOUBLE, rank, BASE+5, MPI_COMM_WORLD);
	  MPI_Send(data->w   , data->N, MPI_DOUBLE, rank, BASE+6, MPI_COMM_WORLD);
	}
	
      }else{ 
	/* master can't send data to itself, so we're using a little trick here.
	 * It first recovers the pointers from the first call ("Parent"), then
	 * shifts the data back to the start, and finally reallocates the array 
	 * to free the unsused memory (i.e. data sent to other cpus). 
	 * Note: memmove() is safer than memcpy().
	 */
	
	dataParent->x = (double *)memmove(dataParent->x, data->x, data->N*NDIM*sizeof(double));
	if(para.corr == GGLENS){
	  dataParent->zerr = (double *)memmove(dataParent->zerr, data->zerr, data->N*sizeof(double));
	  dataParent->e1   = (double *)memmove(dataParent->e1,   data->e1,   data->N*sizeof(double));
	  dataParent->e2   = (double *)memmove(dataParent->e2,   data->e2,   data->N*sizeof(double));
	  dataParent->w    = (double *)memmove(dataParent->w,    data->w ,   data->N*sizeof(double));
	}
	dataParent->x = (double *)realloc(dataParent->x, data->N*NDIM*sizeof(double));	
	if(para.corr == GGLENS){
	  dataParent->zerr = (double *)realloc(dataParent->zerr, data->N*sizeof(double));
	  dataParent->e1   = (double *)realloc(dataParent->e1,   data->N*sizeof(double));
	  dataParent->e2   = (double *)realloc(dataParent->e2,   data->N*sizeof(double));
	  dataParent->w    = (double *)realloc(dataParent->w,    data->N*sizeof(double));
	}
	
	/* copy the new number of points and the dim 
	 * along which they are sorted 
	 */
	dataParent->N   = data->N;
	dataParent->dim = dim;
      }
      
      rank--;
      /* call recursively to transmit the entire data set to all cpus without further splitting */
      if(Ncpu == 0 && rank >= 0) comData(para, data, 0, dim, 0);
    }
  
  }else{
    MPI_Status status;

    MPI_Recv(&data->N, 1, MPI_LONG, MASTER, BASE+0, MPI_COMM_WORLD, &status);
    MPI_Recv(&dim, 1, MPI_INT, MASTER, BASE+1, MPI_COMM_WORLD, &status);
    
    data->x = (double *)malloc(data->N*NDIM*sizeof(double));
    MPI_Recv(data->x, data->N*NDIM, MPI_DOUBLE, MASTER, BASE+2, MPI_COMM_WORLD, &status);
    
    if(para.corr == GGLENS){
      data->zerr = (double *)malloc(data->N*sizeof(double));
      data->e1   = (double *)malloc(data->N*sizeof(double));
      data->e2   = (double *)malloc(data->N*sizeof(double));
      data->w    = (double *)malloc(data->N*sizeof(double));
      MPI_Recv(data->zerr, data->N, MPI_DOUBLE, MASTER, BASE+3, MPI_COMM_WORLD, &status);
      MPI_Recv(data->e1  , data->N, MPI_DOUBLE, MASTER, BASE+4, MPI_COMM_WORLD, &status);
      MPI_Recv(data->e2  , data->N, MPI_DOUBLE, MASTER, BASE+5, MPI_COMM_WORLD, &status);
      MPI_Recv(data->w   , data->N, MPI_DOUBLE, MASTER, BASE+6, MPI_COMM_WORLD, &status);
    }
  }
  
  return;
}
#undef BASE

#define BASE 300
void comResult(const Config para, Result result, long Ncpu, int split){
  /* Sends results back to master if slave. Sums up everything if master. If 
   * split is larger than 0, it means the data has been split before building
   * the tree, hence it requires to sum up the total number of objects as well,
   * but only N2, because N1 is never partitioned. 
   */
  long i, rank;
  
  if(para.rank != MASTER){
    switch(para.corr){
    case AUTO: case CROSS:
      MPI_Send(result.NN, para.nbins*(para.nsamples+1), MPI_LONG, MASTER, BASE+0, MPI_COMM_WORLD);
      MPI_Send(result.N1, para.nsamples+1, MPI_LONG, MASTER, BASE+1, MPI_COMM_WORLD);
      MPI_Send(result.N2, para.nsamples+1, MPI_LONG, MASTER, BASE+2, MPI_COMM_WORLD);
      break;
    case GGLENS:
      MPI_Send(result.GG, para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+0, MPI_COMM_WORLD);
      MPI_Send(result.w,  para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+1, MPI_COMM_WORLD);
      MPI_Send(result.Nsources, para.nbins, MPI_DOUBLE, MASTER, BASE+2, MPI_COMM_WORLD);
      MPI_Send(result.meanR,    para.nbins, MPI_DOUBLE, MASTER, BASE+3, MPI_COMM_WORLD);
      MPI_Send(result.e2,       para.nbins, MPI_DOUBLE, MASTER, BASE+4, MPI_COMM_WORLD);
      break;
    }
  }else{
    
    Result slave;
    MPI_Status status;
      
      switch(para.corr){
      case AUTO: case CROSS:
	slave.NN = (long *)malloc(para.nbins*(para.nsamples+1)*sizeof(long));
	slave.N1 = (long *)malloc((para.nsamples+1)*sizeof(long)); 
	slave.N2 = (long *)malloc((para.nsamples+1)*sizeof(long)); 
	for(rank=1;rank<Ncpu;rank++){
	  MPI_Recv(slave.NN, para.nbins*(para.nsamples+1), MPI_LONG, rank, BASE+0, MPI_COMM_WORLD, &status);
	  MPI_Recv(slave.N1, para.nsamples+1, MPI_LONG, rank, BASE+1, MPI_COMM_WORLD, &status);
	  MPI_Recv(slave.N2, para.nsamples+1, MPI_LONG, rank, BASE+2, MPI_COMM_WORLD, &status);
	  for(i=0;i<para.nbins*(para.nsamples+1);i++) result.NN[i] += slave.NN[i];
	  if(split){
	    /* only N2 may be partitioned */
	    for(i=0;i<para.nsamples+1;i++) result.N2[i] += slave.N2[i];
	  }
	}	
	break;
      case GGLENS:
	slave.GG       = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));
	slave.w        = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));
	slave.Nsources = (double *)malloc(para.nbins*sizeof(double));
	slave.meanR    = (double *)malloc(para.nbins*sizeof(double));
	slave.e2       = (double *)malloc(para.nbins*sizeof(double));
	for(rank=1;rank<Ncpu;rank++){
	  MPI_Recv(slave.GG, para.nbins*(para.nsamples+1), MPI_DOUBLE, rank, BASE+0, MPI_COMM_WORLD, &status);
	  MPI_Recv(slave.w,  para.nbins*(para.nsamples+1), MPI_DOUBLE, rank, BASE+1, MPI_COMM_WORLD, &status);
	  MPI_Recv(slave.Nsources, para.nbins, MPI_DOUBLE, rank, BASE+2, MPI_COMM_WORLD, &status);
	  MPI_Recv(slave.meanR,    para.nbins, MPI_DOUBLE, rank, BASE+3, MPI_COMM_WORLD, &status);
	  MPI_Recv(slave.e2,       para.nbins, MPI_DOUBLE, rank, BASE+4, MPI_COMM_WORLD, &status);
	  for(i=0;i<para.nbins*(para.nsamples+1);i++){
	    result.GG[i] += slave.GG[i];
	    result.w[i]  += slave.w[i];
	  }
	  for(i=0;i<para.nbins;i++){
	    result.Nsources[i] += slave.Nsources[i];
	    result.meanR[i]    += slave.meanR[i];
	    result.e2[i]       += slave.e2[i];
	  }
	}
	break;
      freeResult(para, slave);
    }
  }
}
#undef BASE

/*----------------------------------------------------------------*
 *Configuration                                                   *
 *----------------------------------------------------------------*/

void initPara(int argc, char **argv, Config *para){
  /* Reads the input parameters given in the configuration 
     file and/or in the command line (highest priority) and 
     returns para, containing all the configuration parameters.
     See moredetails in the structure definition in the header file.
  */
  int i, j, noconfigFile = 1;
  long Ncol;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR], list[NFIELD*NCHAR];
  
  /* initialisation */
  srand((unsigned int)time(NULL));
  EPS = determineMachineEpsilon();
  strcpy(MYNAME,"swot");
  
  /* default parameters below */
  strcpy(para->fileInName1,"data1.cat");
  strcpy(para->fileInName2,"data2.cat");
  strcpy(para->fileRanName1,"ran1.cat");
  strcpy(para->fileRanName2,"ran2.cat");
  /* colum ids, default: 
     ids           1  2   3   4   5  6  7 
     [lensing]     RA DEC z  deltaz e1 e2 weight 
     [w(theta)]    RA DEC             
     [wp(rp)]      RA DEC z  deltaz             
     [xi(r)]       X  Y   Z (Mpc) */
  for(i=0;i<NIDSMAX;i++){ 
    para->data1Id[i] = i+1;
    para->data2Id[i] = i+1;
    para->ran1Id[i]  = i+1;
    para->ran2Id[i]  = i+1;
  }
  NDIM            = 2; /* number of dimensions, 2 or 3 */
  para->coordType = RADEC;
  para->distAng   = &distAngSpher;
  para->proj      = THETA;
  para->corr      = AUTO;
  para->estimator = LS;
  para->min       = 0.0001;
  para->max       = 1.0;
  para->nbins     = 20;
  para->log       = 1;
  para->OA        = 0.05;
  para->nsamples  = 32;
  para->err       = JACKKNIFE;
  para->cov_mat   = 0;
  para->deltaz    = 0.03;
  strcpy(para->fileOutName,   "corr.out");
  strcpy(para->fileCovOutName,"cov.out");
  
  /* only master talks */
  if(para->rank == MASTER){
    para->verbose = 1;    
  }else{
    para->verbose = 0;
  }
  
  /* default cosmology <=> WMAP5 */
  para->a[0] = H0;
  para->a[1] = Omega_M;
  para->a[2] = Omega_L;
  para->a[3] = c;
  
  /* ----------------------------------------------------------------------
   *STEP 1: first loop over arguments. Display help or dump default config file */
  for(i=0;i<argc;i++){
    /* help */
    if((!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || argc == 1)){
      if(para->verbose){
      fprintf(stderr,"\n\n\
                          S W O T\n\n\
                (Super W Of Theta) MPI version 0.12\n\n	\
Program to compute two-point correlation functions.\n\
Usage:  %s -c configFile [options]: run the program\n\
        %s -d: display a default configuration file\n\
Important: if using \"RADEC\" coordinates, the angle\n\
in the input catalogues must be in decimal degrees.\n",MYNAME,MYNAME);
      }
      exit(EXIT_FAILURE);
    }
    /* dump a default configuration file ---------------------------- */
    if(!strcmp(argv[i],"-d")){
      printf("#Configuration file for %s\n",MYNAME);
      printf("#----------------------------------------------------------#\n");
      printf("#Input catalogues                                          #\n");
      printf("#----------------------------------------------------------#\n");
      printf("data1          %s  #input data catalogue #1\n",    para->fileInName1);
      printf("cols1          %d,%d\t  #Column ids for data1\n",  para->data1Id[0],para->data1Id[1]);
      printf("data2          %s  #input data catalogue #2\n",    para->fileInName2);
      printf("cols2          %d,%d\t  #Column ids for data2\n",  para->data2Id[0],para->data2Id[1]);
      printf("ran1           %s\t  #input random catalogue #1\n",para->fileRanName1);
      printf("rancols1       %d,%d\t  #Column ids for ran1\n",   para->ran1Id[0],para->ran1Id[1]);
      printf("ran2           %s\t  #input random catalogue #2\n",para->fileRanName2);
      printf("rancols2       %d,%d\t  #Column ids for ran2\n",   para->ran2Id[0],para->ran2Id[1]);
      printf("coord          RADEC\t  #Coordinates: [RADEC,CART,CART3D,RADEC_Z]\n");
      printf("                    \t  #(in degrees if RADEC)\n");
      printf("#----------------------------------------------------------#\n");
      printf("#Correlation options                                       #\n");
      printf("#----------------------------------------------------------#\n");
      printf("corr           auto\t #Type of correlation: [auto,cross,gglens]\n");
      printf("est            ls\t #Estimator [ls,nat,ham]\n");
      printf("range          %g,%g\t #Correlation range. Unity same as \"coord\":\n",para->min,para->max);
      printf("                    \t #in degrees for RADEC, in Mpc otherwise)\n");
      printf("nbins          %d\t #Number of bins\n",para->nbins);
      printf("log            yes\t #Logarithmic bins [yes,no]\n");
      printf("err            jackknife #Resampling method [bootstrap,jackknife]\n");
      printf("nsamples       %d\t #Number of samples for resampling (power of 2)\n",para->nsamples);
      printf("OA             %g\t #Open angle for approximation (value or \"no\") \n",para->OA);
      printf("#----------------------------------------------------------#\n");
      printf("#Cosmology (for gal-gal correlations, w(R) and xi(rp,PI))  #\n");
      printf("#----------------------------------------------------------#\n");
      printf("H0             %g\t #Hubble parameter\n",para->a[0]);
      printf("Omega_M        %g\t #Relative matter density\n",para->a[1]);
      printf("Omega_L        %g\t #Relative energy density (Lambda)\n",para->a[2]);
      printf("deltaz         %g\t #For gg lensing: Zsource > Zlens + deltaz\n",para->deltaz);
      printf("#----------------------------------------------------------#\n");
      printf("#Output options                                            #\n");
      printf("#----------------------------------------------------------#\n");
      printf("proj           theta\t #Axis projection:\n");
      printf("                    \t #if coord = RADEC, [phys,theta]\n");
      printf("                    \t #if coord = CART3D, RADEC_Z [phys] \n");
      printf("                    \t #if corr  = gglens, [como,phys]\n");
      printf("out            %s\t #Output file\n",para->fileOutName);
      printf("cov            no\t #Covariance matrix of errors [yes,no]\n");
      printf("cov_out        %s\t #Covariance output file\n",para->fileCovOutName);
      exit(EXIT_FAILURE);
    }
  }
  
  /* ----------------------------------------------------------------------
   *STEP 2: second loop over arguments. Find config file and read it */
  for(i=0;i<argc;i++){
    if(!strcmp(argv[i],"-c")){
      checkArg(argv[i],argv[i+1],para);
      FILE *configFile = fopenAndCheck(argv[i+1],"r",para->verbose);
      noconfigFile     = 0;
      while(fgets(line,NFIELD*NCHAR,configFile) != NULL){
	if(getStrings(line,item," ",&Ncol)) 
	  setPara(getCharValue(item,1),getCharValue(item,2),para);
      }
      fclose(configFile);
    }
  }
  if(noconfigFile){
    if(para->verbose) fprintf(stderr,"\n%s: **ERROR** config file not found. Exiting...\n", MYNAME);
    exit(EXIT_FAILURE);
  }

  /* ----------------------------------------------------------------------
   *STEP 3: third loop over arguments. Read the command line options                                 
   * (overwrite config file option). */
  for(i=0;i<argc;i++){
    if(*argv[i] == '-') setPara(argv[i]+1,argv[i+1],para);
  }
  
  /* ----------------------------------------------------------------------
   *STEP 4: readjust parameters if needed */
  if(para->corr == GGLENS){
    NDIM  = 3;
    if(para->proj == THETA){
      para->proj = COMO; /*Default projection for gg lensing*/
    }
  }
  
  /* set the angular distance measurement method */
  switch(para->coordType)
    {
    case RADEC:
      para->distAng = &distAngSpher;
      break;
    case RADEC_Z:
      NDIM = 3;
      para->distAng = &distAngSpher;
      break;
    case CART:
      para->distAng = &distAngCart;
      break;
    case CART3D:
      NDIM = 3;
      para->distAng = &dist3D;
      break;
    }
  
  if(para->log){
    para->min = log(para->min);
    para->max = log(para->max);
  }  
  para->Delta = (para->max - para->min)/(double)para->nbins;
   
  return;
}

void setPara(char *field, char *arg, Config *para){
  long j, Ncol;
  char list[NFIELD*NCHAR];
  /* TO DO: write a function with variable number of arguments 
   * see http://en.wikipedia.org/wiki/Stdarg.h 
   * something like:
   * - setArg(para->item1,string,arg)
   * - setArg(para->item2,value,arg)
   * - setArg(para->item3,list,arg,option1,option2,option3)
   * returns 1 if set or 0 if not such as it can be written
   * }else if(setArg(para->item1,string,arg)){}
   * }else if(setArg(para->item1,string,arg)){}
   * }else if(setArg(para->item1,string,arg)){}
   * or
   * set += setArg(para->item1,string,arg)
   * set += setArg(para->item1,string,arg)
   * set += setArg(para->item1,string,arg)
   */
  
  if(!strcmp(field,"data1")){ 
    checkArg(field,arg,para);
    strcpy(para->fileInName1,arg);
  }else if(!strcmp(field,"cols1")){ 
    checkArg(field,arg,para);
    getStrings(arg,list,",",&Ncol);
    for(j=0;j<Ncol;j++)  para->data1Id[j] = getIntValue(list,j+1);
  }else if(!strcmp(field,"data2")){
    checkArg(field,arg,para);
    strcpy(para->fileInName2,arg);
  }else if(!strcmp(field,"cols2")){ 
    checkArg(field,arg,para);
    getStrings(arg,list,",",&Ncol);
    for(j=0;j<Ncol;j++)  para->data2Id[j] = getIntValue(list,j+1);
  }else if(!strcmp(field,"ran1")){
    checkArg(field,arg,para);
    strcpy(para->fileRanName1,arg);
  }else if(!strcmp(field,"rancols1")){ 
    checkArg(field,arg,para);
    getStrings(arg,list,",",&Ncol);
    for(j=0;j<Ncol;j++)  para->ran1Id[j] = getIntValue(list,j+1);
  }else if(!strcmp(field,"ran2")){
    checkArg(field,arg,para);
    strcpy(para->fileRanName2,arg);
  }else if(!strcmp(field,"rancols2")){ 
    checkArg(field,arg,para);
    getStrings(arg,list,",",&Ncol);
    for(j=0;j<Ncol;j++)  para->ran2Id[j] = getIntValue(list,j+1);
  }else if(!strcmp(field,"coord")) {
    checkArg(field,arg,para);
    if(!strcmp(arg,"RADEC"))        para->coordType = RADEC;
    else if(!strcmp(arg,"RADEC_Z"))    para->coordType = RADEC_Z;
    else if(!strcmp(arg,"CART"))    para->coordType = CART;
    else if(!strcmp(arg,"CART3D"))  para->coordType = CART3D;
    else checkArg(field, NULL, para);
  }else if(!strcmp(field,"corr")) {
    checkArg(field, arg, para);
    if(!strcmp(arg,"cross"))       para->corr = CROSS;
    else if(!strcmp(arg,"auto"))   para->corr = AUTO;
    else if(!strcmp(arg,"gglens")) para->corr = GGLENS;
    else checkArg(field, NULL, para);
  }else if(!strcmp(field,"proj")) {
    checkArg(field,arg,para);
    if(!strcmp(arg,"theta")) para->proj = THETA;
    else if(!strcmp(arg,"como"))  para->proj = COMO;
    else if(!strcmp(arg,"phys"))  para->proj = PHYS;
    else checkArg(field, NULL, para);
  }else if(!strcmp(field,"est")) {
    if(!strcmp(arg,"ls"))       para->estimator = LS;
    else if(!strcmp(arg,"nat")) para->estimator = NAT;
    else if(!strcmp(arg,"ham")) para->estimator = HAM;
    else checkArg(field, NULL, para);
  }else if(!strcmp(field,"range")){
    checkArg(field,arg,para);
    getStrings(arg,list,",",&Ncol);
    para->min = getDoubleValue(list,1);
    para->max = getDoubleValue(list,2);
  }else if(!strcmp(field,"nbins")){ 
    checkArg(field,arg,para); 
    para->nbins = atoi(arg);
  }else if(!strcmp(field,"log")){
    checkArg(field,arg,para);
    if(!strcmp(arg,"yes")) para->log = 1;
    else para->log = 0;
  }else if(!strcmp(field,"err")) {
    checkArg(field,arg,para);
    if(!strcmp(arg,"jackknife"))      para->err = JACKKNIFE;
    else if(!strcmp(arg,"bootstrap")) para->err = BOOTSTRAP;
    else checkArg(field, NULL, para);
  }else if(!strcmp(field,"nsamples")){
    checkArg(field,arg,para);
    para->nsamples = atoi(arg);
  }else if(!strcmp(field,"OA")){
    checkArg(field,arg,para);
    if(!strcmp(arg,"no")) para->OA = -1.0;
    else para->OA     = atof(arg);
  }else if(!strcmp(field,"H0")){
    checkArg(field,arg,para);
    para->a[0]   = atof(arg);
  }else if(!strcmp(field,"Omega_M")){
    checkArg(field,arg,para);
    para->a[1]   = atof(arg);
  }else if(!strcmp(field,"Omega_L")){
    checkArg(field,arg,para);
    para->a[2]   = atof(arg);
  }else if(!strcmp(field,"deltaz")){
    checkArg(field,arg,para);
    para->deltaz  = atof(arg);
  }else if(!strcmp(field,"out") || !strcmp(field,"o")){
    checkArg(field,arg,para);
    strcpy(para->fileOutName,arg);
  }else if(!strcmp(field,"cov")){
    checkArg(field,arg,para);
    if(!strcmp(arg,"yes")) para->cov_mat = 1;
    else para->cov_mat = 0;
  }else if(!strcmp(field,"cov_out")){
    checkArg(field,arg,para);
    strcpy(para->fileCovOutName,arg);
  }else if(!strcmp(field,"c")){
    /* do nothing, this is the option for the config file */
  }else{
    if(para->verbose) fprintf(stderr,"%s: **ERROR** %s is not a valid option. Exiting...\n", MYNAME, field);
    
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  
  return;
}

void checkArg(char *field, char *arg, Config *para){
  /* Checks if the option has a valid argument. */
  if(arg == NULL || *arg == '-'){
    if(para->verbose) fprintf(stderr,"%s: **ERROR** option %s has no valid argument. Exiting...\n", MYNAME, field);
    
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  return;
}

/*----------------------------------------------------------------*
 *Utils - Numeric                                                 *
 *----------------------------------------------------------------*/

double distComo(double z, const double a[4]){
  /* Return the comoving distance of z, given the cosmological 
   * parameters a. */
  
  int n = 1000;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc (n);
  
  gsl_function F;
  F.function = &drdz;
  F.params   = (void *)a;

  double result, error;
  gsl_integration_qag(&F, 0.0, z, 0, 1e-7, n, 6, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return result;
}

double drdz(double x, void * params){
  /* Inverse of E(z). With following cosmological parameters:
   * a[0] = H0;
   * a[1] = Omega_M;
   * a[2] = Omega_L;
   * a[3] = c;
   */
  
  double *a = (double *) params;
  return a[3]/(a[0]*sqrt(a[1]*pow(1+x,3.0)+a[2]));
}

void splitData(const Config para, const Point data, int dim, Point *dataLeft, Point *dataRight){
  /* Takes a full array of data and sorts it in place along 
   * the "dim" coordinate. Then returns the address of the two 
   * children in dataLeft and dataRight, and the number of 
   * points in NLeft and NRight. Since the partition is 
   * made "in place" it does not use any extra memory.
   */
  
  /* sort data along the "dim" coordinate */
  quickSort(para, data, dim, 0, data.N - 1);
  
  /* if N is odd, puts more objects on the left (pure convention) */
  switch(PARITY(data.N)){ 
  case EVEN: dataLeft->N = data.N/2;     dataRight->N = data.N/2;     break;
  case ODD:  dataLeft->N = (data.N+1)/2; dataRight->N = (data.N-1)/2; break;
  }
  
  dataLeft->dim = dim ; dataRight->dim = dim;           /* sorting dimension                    */
  copyPointAddress(para, dataLeft,  data, 0);           /* left child.  starts at data + 0      */
  copyPointAddress(para, dataRight, data, dataLeft->N); /* right child. starts at data + NLeft  */
  
  return;
}


void quickSort(const Config para, Point data, int dim, long start, long end){
  /* Sorts data in place (increasing order) along the dimension "dim". 
   * See http://en.wikipedia.org/wiki/Quicksort and Numerical Recipes, Press et al.
   * There's room for imporvement here. See qsort.c in particular. Also maybe 
   * think of replacing functions by macros.
   * ATTENTION: when coded recursively (albeit more beautiful), 
   * it's supposed to be a bit slower (Jean: I haven't tested it).
   */
  
  long i, j, middle, jstack=-1, Nmin = 7, stack[NSTACK];
  double pivotValue;
  Point d = createPoint(para,1);
  
  /* first set the dimension along which points are going to be sorted */
  data.dim = dim;

  /* sanity check */
  if(end - start + 1 > data.N){
    if(para.verbose) fprintf(stderr,"%s: **ERROR** range width (%zd) in quickSort() larger than the number of points (%zd). Exiting...\n", MYNAME, end - start + 1, data.N);
    exit(EXIT_FAILURE);
  }
  
  for(;;){
    if(end - start + 1 <= Nmin){    
      /* perform insertion sort if the size of the array is small (<= Nmin).
       * Nmin = 7 is usually recommended but word on the street says
       * it depends on the platform (it's not very important here). */
      
      for (j=start+1;j<=end;j++) {
	copyPoint(para, d, 0, data, j);
      	pivotValue = data.x[NDIM*j+dim];
      	for (i=j-1;i>=start;i--) {
      	  if (data.x[NDIM*i+dim] < pivotValue ) break;
	  copyPoint(para, data, i+1, data, i);
      	}
	copyPoint(para, data, i+1, d, 0);
      }
      
      /* if the stack is empty, the loop has visited all
       * the data table. Job over. */
      if(jstack < 0) break;
      
      /* next turn */
      end   = stack[jstack--];
      start = stack[jstack--];
      
    }else{
      /* qsort algorithm (see NR) */
      middle=(start+end)/2;
      swapPoint(para, data, start+1, middle);

      /* choose median of left, center, and right element 
       * as partitioning element d. Also rearrange so that
       * data[start] < data[start+1] < data[end] */
      if(data.x[NDIM*start+dim] > data.x[NDIM*end+dim])
	swapPoint(para, data, start, end);
      if(data.x[NDIM*(start+1)+dim] > data.x[NDIM*end+dim]) 
	swapPoint(para, data, start+1, end);
      if(data.x[NDIM*(start)+dim] > data.x[NDIM*(start+1)+dim]) 
	swapPoint(para, data, start, start+1);      

      /* initialize pointers for partitioning */
      i = start+1;
      j = end;

      pivotValue = data.x[NDIM*(start+1)+dim];
      copyPoint(para, d, 0, data, start+1);
      
      /* beginning of innermost loop */
      for(;;) {
	do i++;	while (data.x[NDIM*i+dim] < pivotValue);
	do j--; while (data.x[NDIM*j+dim] > pivotValue);
	if (j < i) break;
	swapPoint(para, data, i, j);
      }
      /* end of innermost loop */
      
      /* insert partitioning element */
      copyPoint(para, data, start+1, data, j);
      copyPoint(para, data, j, d, 0);
      /* increment stack */
      jstack += 2;

      /* if more than 2^NSTACK points */
      if(jstack >= NSTACK){
	if(para.verbose) fprintf(stderr,"%s: **ERROR** NSTACK (%d) too small in quickSort. Exiting...\n", MYNAME, NSTACK);
	exit(EXIT_FAILURE);
      }
      
      /* for next run */
      if(end - i+1 >= j-start){
	stack[jstack]   = end;
	stack[jstack-1] = i;
	end = j-1;
      }else{
	stack[jstack]   = j-1;
	stack[jstack-1] = start;
	start = i;
      }      
    }
  }
  
  freePoint(para, d);
  return;
}

void copyPoint(const Config para, Point a, long i, Point b, long j){
  /* Copies values from b into a. a and b can be the same object.
   * CAUTION: the destination "a" might not be sorted nor
   * filled with the same number of points. That's why it is 
   * only a pair-wise copy (and not Point.N or Point.dim). */
  long dim;
  
  for(dim=0;dim<NDIM;dim++) a.x[NDIM*i+dim] = b.x[NDIM*j+dim];
  if(para.corr == GGLENS){
    a.zerr[i] = b.zerr[j];
    a.e1[i]   = b.e1[j];
    a.e2[i]   = b.e2[j];
    a.w[i]    = b.w[j];
  }
  
  return;
}


void getMeanPoint(const Config para, Point a, long i, Point point){
  /* get the mean properties of an ensemble of 
   * points "point" and put it in a(i) */
  
  /* TO DO: for zerr, it should be the dispersion of point 
   * and not the mean error. */
  
  
  long n, dim;
  /* set to zero */
  for(dim=0;dim<NDIM;dim++) a.x[NDIM*i+dim] = 0.0;
  if(para.corr == GGLENS){
    a.zerr[i] = 0.0;
    a.e1[i]   = 0.0;
    a.e2[i]   = 0.0;
    a.w[i]    = 0.0;
  }
  /* sum up */
  for(n=0;n<point.N;n++){
    for(dim=0;dim<NDIM;dim++) a.x[NDIM*i+dim] += point.x[NDIM*n+dim];
    if(para.corr == GGLENS){
      a.zerr[i] +=  point.zerr[n];
      a.e1[i]   +=  point.e1[n];
      a.e2[i]   +=  point.e2[n];
      a.w[i]    +=  point.w[n];
    }
  }
  /* divide by N */
  for(dim=0;dim<NDIM;dim++) a.x[NDIM*i+dim] /= (double)point.N;
  if(para.corr == GGLENS){
    a.zerr[i] /= (double)point.N;
    a.e1[i]   /= (double)point.N;
    a.e2[i]   /= (double)point.N;
    a.w[i]    /= (double)point.N;
  }
  
  return;
}

void copyPointAddress(const Config para, Point *a, const Point b, long shift){
  /* Copies the array addresses of b into a. ATTENTION: as it does not 
   * copy the value but the address, it means that any change to a will
   * also affect b. Set shift to 0 to start at the beginning of b or 
   * any number to access points at "shift", e.g. b.*[shift]. 
   *
   * Not familiar with pointers? Read this: http://www.cs.cf.ac.uk/Dave/C/node10.html 
   */
  long dim;
  
  a->x = b.x + shift*NDIM;
  if(para.corr == GGLENS){
    a->zerr = b.zerr + shift;
    a->e1   = b.e1   + shift;
    a->e2   = b.e2   + shift;
    a->w    = b.w    + shift;
  }
  
  return;
}


void swapPoint(const Config para, Point point,long i, long j){
  /* Swaps i and j in the array "point". */

  long dim;
  double tmp;
  
  for(dim=0;dim<NDIM;dim++){
    tmp                 = point.x[NDIM*i+dim];
    point.x[NDIM*i+dim] = point.x[NDIM*j+dim];
    point.x[NDIM*j+dim] = tmp;
  }
  if(para.corr == GGLENS){
    tmp           = point.zerr[i];
    point.zerr[i] = point.zerr[j];
    point.zerr[j] = tmp;
    tmp           = point.e1[i];
    point.e1[i]   = point.e1[j];
    point.e1[j]   = tmp;
    tmp           = point.e2[i];
    point.e2[i]   = point.e2[j];
    point.e2[j]   = tmp;
    tmp           = point.w[i];
    point.w[i]    = point.w[j];
    point.w[j]    = tmp;
  }
  
  return;
}


Point createPoint(const Config para, long N){
  Point point;

  point.N   = N;
  point.dim = -1;  
  
  point.x = (double *)malloc(N*NDIM*sizeof(double));
  if(para.corr == GGLENS){
    point.zerr = (double *)malloc(N*sizeof(double));
    point.e1   = (double *)malloc(N*sizeof(double));
    point.e2   = (double *)malloc(N*sizeof(double));
    point.w    = (double *)malloc(N*sizeof(double));
  }
  
  return point;
}

void freePoint(const Config para, Point point){
  
  free(point.x);
  if(para.corr == GGLENS){
    free(point.zerr);
    free(point.e1);
    free(point.e2);
    free(point.w);
  }
  
  return;
}



Point readCat(const Config para, char *fileInName, int id[NIDSMAX]){
  /* Reads fileIn and returns the data array 
   * (see header file for details about the 
   * data structure.
   */
  long n, N, dim, Ncol;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR];
  
  FILE *fileIn = fopenAndCheck(fileInName, "r", para.verbose);
  
  /* get size of file and allocate data */
  N = 0;
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL)
    if(*(line) != '#' && *(line) != '\0' && *(line) != '\n') N++;
  rewind(fileIn);
  
  
  if(N < para.size){
    /* if no point found in fileIn or if the number of 
     * points is less than the number of cpu (who knows...) 
     */
    if(para.verbose) fprintf(stderr,"\n%s: **ERROR** no valid point found in the file or too many cpus. Exiting...\n", MYNAME);
    //  exit(EXIT_FAILURE);      
  }
  
  Point data = createPoint(para, N);
  
  n = 0;
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){
    
    if(getStrings(line,item," ",&Ncol)){
      /* read the value at the column definied by "id[i]" */
      for(dim=0;dim<NDIM;dim++) data.x[NDIM*n+dim] = getDoubleValue(item,id[dim]);
      if(para.corr == GGLENS){
	Ncol -= 4;
	//Check for positive values here ??
	//TO DO add deltaz if definied in para file
	data.zerr[n] = getDoubleValue(item,id[NDIM+0]);
	data.e1[n]   = getDoubleValue(item,id[NDIM+1]);
	data.e2[n]   = getDoubleValue(item,id[NDIM+2]);
	data.w[n]    = getDoubleValue(item,id[NDIM+3]);
      }
      if(Ncol < NDIM){
	/*
	  if(para.verbose){
	  fprintf(stderr,"\n%s: **ERROR** (row %zd), the number of columns is less than expected. Check configuration or input file. Exiting...\n", MYNAME, n+1);
	  }
	
	exit(EXIT_FAILURE);      
	*/
      }      
      /* increment n */
      n++;
    }
  }
  
  fclose(fileIn);
  return data;
}

FILE *fopenAndCheck(const char *fileName, char *mode, int verbose){
  /* Checks if fileName exists and opens it. Exits otherwise. */
    FILE *fileTmp = fopen(fileName, mode);
    if (fileTmp == NULL){
      if(verbose) printf("%s: **ERROR** %s not found. Exiting...\n",MYNAME,fileName);
      //TO DO proper error display
      exit(EXIT_FAILURE); 
  }
    return fileTmp;
}

int getStrings(char *line, char *strings, char *delimit, long *N){
  /* Extract each word/number in line separated by delimit 
   * and returns the array of items in strings.
   */
  int i,j,begin,length;
  
  if(line == NULL || line[0] == '\n' || line[0] == '#' || line[0] == '\0') return 0;
  
  i = 0;
  j = 0;
  while(line[i] != '\0' && line[i] != '#' && line[i] != '\n'){
    begin = i;
    while(line[i] == *delimit || line[i] == '\t' && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
    begin = i;
    while(line[i] != *delimit && line[i] != '\t' && line[i] != '\0' && line[i] != '#' && line[i] != '\n') i++;
    length = i - begin;
    if(length > 0){
      strncpy(strings+NCHAR*j,&line[begin],length);
      strcpy(strings+NCHAR*j+length,"\0");
      j++;
    }
  }
 
  (*N) = j;

  if(*N > 0){
    return 1;
  }else{
    return 0;
  }
}

double determineMachineEpsilon(){
  double u, den;
  u = 1.0;
  do {
    u /= 2.0;
    den = 1.0 + u;
  } while(den>1.0);
  
  return(10.0 * u);
}

long randInt(long N){
  /* Returns an random integer between 0 and N-1. */
  return (long)((double)rand()/((double)RAND_MAX + 1)*N);
}

void printCount(const long count, const long total, const long step, int verbose){
  if(verbose){
    if(count%step == 0){
      fflush(stdout);
      fprintf(stderr,"\b\b\b\b\b\b\b%6.2f%%",100.0*(double)count/(double)total);
    }
  }
}

void comment(const Config para, char *commentString){
  if(para.verbose){
    fflush(stderr); 
    fprintf(stderr,commentString);
  }
  return;
}

long lowestPowerOfTwo(long n, long *pow){
  long i = 1;
  
  *pow = 0;
  
  if(n == 0) return 0;
  
  for(;;){
    if(i == n){
      break;
    }else if(i > n){
      i >>= 1;
      (*pow)--;
      break;
    }
    i <<= 1;
    (*pow)++;
  }
  
  return i;
  
}

/*----------------------------------------------------------------*
 *obsolete                                                        *
 *----------------------------------------------------------------*/

#define leaf(tree,node) ((tree)->left[(node)] == -1)
#define node(tree,node) ((tree)->left[(node)] > -1)

Result Npairs2(const Config *para, const Tree *tree1, const long root1, const Tree *tree2, const long root2){ 
  /* returns the number of pairs if tree1 != tree2, and 2 times 
   * the number of pairs if tree1 == tree2, so that the treatment of pairs 
   * in w(theta) estimators is identical (numerical trick). As a result
   * error bars must be divided by sqrt(2.0).
   */
  
  long N, k, count = 0;
  Result result;
  result.NN = (long *)malloc(para->nbins*sizeof(long)); 
  for(k = 0; k < para->nbins; k++) result.NN[k] = 0;  
  
  double d;
  
  long i, j, stackIndex1 = 0, stackIndex2 = 0, stack1[NSTACK], stack2[NSTACK];
  i = root1;
  for(;;){
    j = root2;
    for(;;){
      /* computes only half of the pairs for auto-correlation */
      if(tree1 == tree2 && i < j) goto skip;
      d = para->distAng(tree1, &i, tree2, &j);
      if(node(tree1,i) && tree1->r[i]/d > para->OA){   /* go down tree1        */
	stack1[stackIndex1++] = tree1->right[i];       /* put right node on the stack1  */
	i                     = tree1->left[i];        /* process left node immediately */
	if(node(tree2,j) && tree2->r[j]/d > para->OA){ /* go down tree2      */
	  stack2[stackIndex2++] = tree2->right[j];
	  j                     = tree2->left[j];
	}
      }else if(node(tree2,j) && tree2->r[j]/d > para->OA){ /* go down tree2  */
	stack2[stackIndex2++] = tree2->right[j];
	j                     = tree2->left[j];
      }else{
	//if(tree1 == tree2 && i == j) goto skip;
	if(para->log) d = log(d);
	k = floor((d - para->min)/para->Delta);
	if(k < para->nbins){
	  if(k >= 0) result.NN[k] += tree1->N[i]*tree2->N[j];
	}
      skip: if(stackIndex2 == 0) break;  /* has visited all nodes in tree2..       */        
	j = stack2[--stackIndex2];       /* ...or will process next on the stack   */       
      }
    }
    if(stackIndex1 == 0) break;
    i = stack1[--stackIndex1];
    printCount(i+1,tree1->size,100,para->verbose);
  }
  
  if(para->verbose) printf("\b\b\b\b\b\b\b%6.2f%%\n",100.0); 
  
  if(tree1 == tree2){
    /* ATTENTION: this is not accurate if the number of objects is small.
     * [Jean] I was thinking that a way to deal with it would be to
     * check whether sqrt(1+8*NN) is rational or not. If it is then
     * NN *= 2 is correct; if not, n = sqrt(2*NN) and NN(true) = n*(n-1)/2. 
     * I'm wondering if this is a problem or not, given that we 
     * deal with large samples. However at small scales we have very few 
     * objects, so it could be a problem.*/
       for(k = 0; k < para->nbins; k++) result.NN[k] *= 2;
  }
  
  return result;  
  
  /* Alternative algorithms:
  


  i = ROOT;
  for(;;){
    j = ROOT;
    for(;;){
      if(node(tree1,i) && node(tree2, j)){       // go down tree1 and tree2  
  	stack1[stackIndex1++] = tree1->right[i];  // put right node on the stack
  	i = tree1->left[i];                       // process left node first    
  	stack2[stackIndex2++] = tree2->right[j];
  	j = tree2->left[j];
      }else if(node(tree1,i)){                   // go down tree1              
  	stack1[stackIndex1++] = tree1->right[i];
  	i = tree1->left[i];
      }else if(node(tree2,j)){                   // go down tree2              
  	stack2[stackIndex2++] = tree2->right[j];
  	j = tree2->left[j];
      }else{                                     // correlate nodes            
  	result.NN[0] += tree1->N[i]*tree2->N[j];
  	if(stackIndex2 == 0) break;              // visited all nodes in tree2->..
       	j = stack2[--stackIndex2];               // ...or process next on the stack 
      }
    }
    if(stackIndex1==0) break;  // visited all nodes in tree1 and tree2-> job over...
    i = stack1[--stackIndex1]; // ...or process next on the stack                  
    printCount(i+1,tree1->size,100,para->verbose);
  }
  
  //when node condition is independent from each other ---------------- 
  for(;;){
    printCount(&i,&tree1->size,100,para->verbose);
    if(leaf(tree1,i)){
      j = 0;
      stackIndex2 = 0;
      for(;;){
  	if(leaf(tree2,j)){
  	  result.NN[0] += tree1->N[i]*tree2->N[j];
  	  stackIndex2--;
  	  if(stackIndex2<0) break;
  	  j = stack2[stackIndex2];
  	}else{
  	  stack2[stackIndex2++] = tree2->right[j];
  	  j = tree2->left[j];
  	}
      }
      stackIndex1--;
      if(stackIndex1<0) break;
      i = stack1[stackIndex1];
    }else{
      stack1[stackIndex1++] = tree1->right[i];
      i = tree1->left[i];
    }
  }    


  */
  
  //if(para->verbose) printf("\b\b\b\b\b\b\b%6.2f%%\n",100.0); 
  
  //return result;
}

