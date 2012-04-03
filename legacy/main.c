/*-------------------------------------------------------------------------*
 *swot (Super W Of Theta)  mpi version                                     *
 *Jean Coupon, Alexie Leauthaud (2012)                                     *
 *-------------------------------------------------------------------------*
 
Program to compute two-point correlation functions. Supports	
auto and cross correlations and galaxy-galaxy lensing.

TO DO:

- reduce memory usage: create "light" node if rmax/dmax = OA, 
  float for coordinates, warning in case of large cats?
- optimize DOL/DOS, can put some of this in memory, etc ...
- Make an option to save binary file of source tree (make an 
  option to read it in next time..). Only really need to compute 
  the source tree once because this never changes.
- Fast option version for randoms: for randoms we don't need e2, 
  nsource, don't need r_moyen ...Also randoms does not need bootstrap.
- For randoms: read in a list of 50 different random files (these will 
  be pre-computed). Only calcualte the source tree once. Only calcualte 
  DOL once (randoms will all have same set of z).
- Do the randoms one by one or all together at once? -> Still need 
  to decide this question.
- For the randoms: will give a **list** of files in input; for each 
  random computed, writes output: [long term requirement if a single sample 
  takes ages to be computed] store nodes completed + sum of pairs  every day ?  
- Long term: bootstrap over regions (instead of Lens by Lens). 
  Jackknife: remove one node at a time. Do RA/DEc boxes versus 
  RA/DEC/z boxes? But do boxes have to be independant? 
- area bootstrap -> weigh nodes instead of objects (e.g. 0,0,1,0). 
  Computed sum of nodes X weights in the end. 
- merge option reading from command line and config file
- options to take into account the East-West orientation
- measure the tangential shear

Contributions:
- the algorithm to compute the number of pairs from a tree is 
  based on, but slightly different to Martin Kilbinger's Ahtena code:
  http://www2.iap.fr/users/kilbinge/athena/
- the gal-gal lensing algorithm is based on
  Alexie Leauthaud's code 
  (see  Leauthaud et al. (2010),  2010ApJ...709...97L).

v 1.42 March 7th 2012 [Jean]
- error on photo-z for sources now
  taken into account

v 1.41 Feb 19th 2012 [Alexie]
- fixed a small bug with recording e2

v 1.4 Feb 14th 2012 [Jean]
- now displays <R>, e2 and Nsources
- output file format is improved

v 1.34 Jan 18th 2012 [Alexie + Jean]
- Added comments to code

v 1.33 Jan 10th 2012 [Jean]
- Fixed a bug to properly sum up the number of 
  sources per bin from all cpus.
- Fixed a bug to get phys or comoving coordinates (Alexie)

v 1.3 Dec 2011 [Alexie]
- Fixed the gg-lensing printf statement to not divide by zero when 
there are no source galaxies.
- Added a printout statement to print the number of sources per bin.

v 1.2 Dec 2011
- added the option -proj to give physical projection 
  for w(R) and Delta sigma

v 1.1 Nov 2011
- added the photo-z error in input catalogue
  lensing signal is added for zs > zl + sigz(zl)+sigz(zs)
  
v 1.0 Nov 2011
- first version, from old "swot" code + parallel 
  computing support (mpi) and gal-gal lensing.
  
*/

#include "main.h"

int main(argc,argv)
     int argc;
     char **argv;
{ 
  int  master = 0, firstCall = 1, verbose = 0;
  int rank, size, nboots;
  
  /* MPI initialization */
  MPI_Status status;
  MPI_Init(&argc, &argv);
  /* "rank" is the id of the current cpu in use [0: master, size-1: last cpu]    */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  /* "size" is the total number of cpus in use. "-np" option when calling mpirun */
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
  
  /* initialization: all cpus read config file */
  Config para; readPara(argc,argv,&para);
  
  /* last rank likely to be the slowest and used to display the messages */
  if(rank == size - 1) verbose = 1;
  double t0 = MPI_Wtime();
  
  switch (para.corr){
  case AUTO:
    /* two-point autocorrelation function */
    autoCorr(para,rank,size,verbose);
    break;
  case CROSS:
    /* two-point cross-correlation function */
    crossCorr(para,rank,size,verbose);
    break;
  case GGLENS:
    ggCorr(para,rank,size,verbose);
    /* galaxy-galaxy lensing two-point cross-correlation function */
    break;
  }
  
  /* computation time */
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == master) printf("Elapsed time: %5.2f s\n", MPI_Wtime() - t0);
  
  MPI_Finalize();
  return EXIT_SUCCESS;
}

/*----------------------------------------------------------------*
 *Main routines                                                   *
 *----------------------------------------------------------------*/


int ggCorr(Config para, int rank, int size, int verbose){
  /*Computes the galaxy-galaxy shear correlation function*/
  long i,j,k,n,Nsources, Nlenses;
  int  nboots,master = 0, firstCall = 1;
  MPI_Status status;
 
  /* tree for lenses */
  Node *lensTree   = readAndCreateTree(para,para.fileInName1,para.data1Id,rank,&Nlenses);
  nboots           = para.nboots; 
  
  /* tree for sources <- this is a memory monster */
  para.nboots      = 0;      /* no bootstrap for sources (errors from weights) */
  Node *sourceTree = readAndCreateTree(para,para.fileInName2,para.data2Id,rank,&Nsources);
  
  para.nboots      = nboots; /* re-set bootstrap for lenses */
  
  /* initialization. para.nboots+1: because first bin is without bootstrap (thus +1) */
  double *SigR          = (double *)malloc(para.nbins*(para.nboots+1)*sizeof(double));
  double *weight        = (double *)malloc(para.nbins*(para.nboots+1)*sizeof(double));
  /* "info" contains the useful pieces of information:
     0 -> number of sources 
     1 -> mean R
     2 -> e2
  */
  para.Ninfo   = 3;
  double *info = (double *)malloc(para.Ninfo*para.nbins*sizeof(double));
  
  for(i=0;i<para.nbins;i++){
    for(n=0;n<para.nboots+1;n++){
      weight[n+(para.nboots+1)*i] = 0.0;
      SigR[n+(para.nboots+1)*i]   = 0.0;
    }
    for(j=0;j<para.Ninfo;j++){
      info[j+(para.Ninfo)*i]      = 0.0;
    }
  } 
  
  /* Check point for processors (all processors wait here)*/
  MPI_Barrier(MPI_COMM_WORLD);
  if (verbose) printf("Correlating lenses with sources...       ");  
  
  /* calling gg-lensing: compute cross-correlation here */
  gg(&para,lensTree,sourceTree,SigR,weight,info,rank,size,firstCall,verbose);
  
  /* send "SigR", "weight" and "info" to the master processor. "1" and "2" here are just an ID number */ 
  /* note: when using MPI_Send or Recieve... easier to send single type tables (not structures)       */
  MPI_Send(SigR,   para.nbins*(para.nboots+1), MPI_DOUBLE, master, 1, MPI_COMM_WORLD); /*1 here is an ID for SigR  */
  MPI_Send(weight, para.nbins*(para.nboots+1), MPI_DOUBLE, master, 2, MPI_COMM_WORLD); /*2 here is an ID for weight*/
  MPI_Send(info,   para.nbins*para.Ninfo,      MPI_DOUBLE, master, 3, MPI_COMM_WORLD); /*3 here is an ID for info  */
  
  /* summing over contributions from all cpus. From now on only the master is working */
  if(rank == master){
    /* the following is only used by master */
    double *SigRRank   = (double *)malloc(para.nbins*(para.nboots+1)*sizeof(double)); 
    double *weightRank = (double *)malloc(para.nbins*(para.nboots+1)*sizeof(double));
    double *infoRank   = (double *)malloc(para.nbins*para.Ninfo*sizeof(double));
    
    /* loop over cpus, "size" is N cpus, start at 1 because 0 is master*/
    for (k = 1; k < size; k++) { 
      /* master receives information from "k" and stores the info in [...]Rank */
      MPI_Recv(SigRRank,   para.nbins*(para.nboots+1), MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status); 
      MPI_Recv(weightRank, para.nbins*(para.nboots+1), MPI_DOUBLE, k, 2, MPI_COMM_WORLD, &status); 
      MPI_Recv(infoRank,   para.nbins*para.Ninfo     , MPI_DOUBLE, k, 3, MPI_COMM_WORLD, &status); 
      /* add together contributions from processors : */
      for(i = 0;i<para.nbins;i++){
	for(n = 0;n<para.nboots+1;n++){
	  SigR[n+(para.nboots+1)*i]   += SigRRank[n+(para.nboots+1)*i];
	  weight[n+(para.nboots+1)*i] += weightRank[n+(para.nboots+1)*i];
	}
	for(j = 0;j<para.Ninfo;j++){
	  info[j+para.Ninfo*i] += infoRank[j+para.Ninfo*i];
	}
      }
    }
    free(SigRRank);
    free(weightRank);
    free(infoRank);
   
    
    /* R(Mpc) */
    /* TO CHANGE: calculte mean here, get sum of R's and divide by N_tot_source later ...
    double *R = (double *)malloc(para.nbins*sizeof(double));
    if(para.log){
      for(i=0;i<para.nbins;i++) R[i] = exp(para.min+para.Delta*(double)i+para.Delta/2.0);  
    } else {
      for(i=0;i<para.nbins;i++) R[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
    }
    */

    /* SigR, only calculate if more than 3 sources (avoid division by zero) */
    for(i=0;i<para.nbins;i++){
      for(n=0;n<para.nboots+1;n++){
	if(info[0+para.Ninfo*i] > 3.0) SigR[n+(para.nboots+1)*i] /= weight[n+(para.nboots+1)*i];
      }
    }
    
    /* bootstrap errors */
    // TO DO: bootstrap on RA/DEC boxes and write covar matrix
    double *SigR_err_boot = (double *)malloc(para.nbins*sizeof(double));
    for(i=0;i<para.nbins;i++){
      SigR_err_boot[i] = 0.0;
      for(n=1;n<para.nboots+1;n++){
	SigR_err_boot[i] += (SigR[(para.nboots+1)*i] - SigR[n+(para.nboots+1)*i])*(SigR[(para.nboots+1)*i] - SigR[n+(para.nboots+1)*i]);
      }
      if(para.nboots > 1) SigR_err_boot[i] = sqrt(SigR_err_boot[i]/(double)(para.nboots-1));
    }
    
    /* print out gg-lensing results */
    double R,Rmean;
    FILE *fileOut = fopen(para.fileOutName,"w");
    fprintf(fileOut,"#gal-gal lensing. Sigma(R), linear approximation\n");
    if(para.proj == COMO){
      fprintf(fileOut,"#R(Mpc) in comoving coordinates\n");
    }else if(para.proj == PHYS){
      fprintf(fileOut,"#R(Mpc) in physical coordinates\n");
    }
    fprintf(fileOut,"#Nbootstraps = %d\n",para.nboots);
    fprintf(fileOut,"#  R(Mpc) SigR(Msun/pc^2)  Err_weights   Err_bootstrap     Nsources       <R>         e2\n");
    for(i=0;i<para.nbins;i++){
      /* don't divide by zero if there are too few objects in the bin */
      if(info[0+para.Ninfo*i] < 3.0) weight[(para.nboots+1)*i] = 1.0;
      /* R and Rmean (weighted) */
      if(para.log){ 
	R     = exp(para.min+para.Delta*(double)i+para.Delta/2.0);
	Rmean = exp(info[1+para.Ninfo*i]/weight[(para.nboots+1)*i]);
      }else{
	R     = para.min+para.Delta*(double)i+para.Delta/2.0;
	Rmean = info[1+para.Ninfo*i]/weight[(para.nboots+1)*i];
      }
      fprintf(fileOut,"%12.7f %12.7f %12.7f %12.7f %15zd %12.7f %12.7f\n", 
	      R, -SigR[(para.nboots+1)*i],				\
	      sqrt(1.0/weight[(para.nboots+1)*i]), SigR_err_boot[i],	\
	      (long)info[0+para.Ninfo*i],				\
	      Rmean,							\
	      info[2+para.Ninfo*i]/weight[(para.nboots+1)*i]);
    }
    fclose(fileOut);
    free(SigR_err_boot);
  }
  
  free(weight);
  free(SigR);
  free(info);
  return SUCCESS;
}

int autoCorr(Config para, int rank, int size, int verbose){
  /* Computes the auto correlation function */
  long i,j,k,Ndata, Nrandom, *DD, *DR, *RR;
  int  master = 0, firstCall = 1;
  MPI_Status status;
  
  Node *dataTree   = readAndCreateTree(para,para.fileInName1,para.data1Id,rank,&Ndata);
  Node *randomTree = readAndCreateTree(para,para.fileRanName1,para.ran1Id,rank,&Nrandom);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(verbose) printf("DD...       ");  DD = Npairs(&para,dataTree,  dataTree,   rank, size, firstCall,verbose);
  if(verbose) printf("DR...       ");  DR = Npairs(&para,dataTree,  randomTree, rank, size, firstCall,verbose); free_Node(dataTree);
  if(verbose) printf("RR...       ");  RR = Npairs(&para,randomTree,randomTree, rank, size, firstCall,verbose); free_Node(randomTree);  

  MPI_Send(DD,  para.nbins*(para.nboots+1), MPI_LONG, master, 1, MPI_COMM_WORLD);
  MPI_Send(DR,  para.nbins*(para.nboots+1), MPI_LONG, master, 2, MPI_COMM_WORLD);
  MPI_Send(RR,  para.nbins*(para.nboots+1), MPI_LONG, master, 3, MPI_COMM_WORLD);
  
  if(rank == master){
    long *DDRank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    long *DRRank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    long *RRRank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    
    for (k = 1; k < size; k++) {
      MPI_Recv(DDRank, para.nbins*(para.nboots+1), MPI_LONG , k, 1, MPI_COMM_WORLD, &status);  
      MPI_Recv(DRRank, para.nbins*(para.nboots+1), MPI_LONG , k, 2, MPI_COMM_WORLD, &status);  
      MPI_Recv(RRRank, para.nbins*(para.nboots+1), MPI_LONG , k, 3, MPI_COMM_WORLD, &status);  
      for(i = 0;i<para.nbins;i++){
	for(j=0;j<para.nboots+1;j++){
	  DD[j+(para.nboots+1)*i] += DDRank[j+(para.nboots+1)*i];
	  DR[j+(para.nboots+1)*i] += DRRank[j+(para.nboots+1)*i];
	  RR[j+(para.nboots+1)*i] += RRRank[j+(para.nboots+1)*i];
	}
      }
    }
    free(DDRank);
    free(DRRank);
    free(RRRank);
  }
  
  if(rank == master)
    switch(para.estimator){
    case LS:
      writeLandySzalay(para,DD,DR,DR,RR,Ndata,Nrandom,Ndata,Nrandom);
      break;
    case NAT:
      writeNatural(para,DD,RR,Ndata,Nrandom,Ndata,Nrandom);
      break;
    case HAM:
      writeHamilton(para,DD,DR,DR,RR,Ndata,Nrandom,Ndata,Nrandom);
      break;
    }

  free(DD);
  free(DR);
  free(RR);
  
  return SUCCESS;
}

int crossCorr(Config para, int rank, int size, int verbose){
  /*Computes the auto correlation function*/
  long i,j,k, Ndata1, Nrandom1, Ndata2, Nrandom2, *D1D2, *D1R1, *D2R2, *R1R2;
  int  master = 0, firstCall = 1;
  MPI_Status status;
  
  Node *data1Tree   = readAndCreateTree(para,para.fileInName1,para.data1Id,rank,&Ndata1);
  Node *data2Tree   = readAndCreateTree(para,para.fileInName2,para.data2Id,rank,&Ndata2);
  Node *random1Tree = readAndCreateTree(para,para.fileRanName1,para.ran1Id,rank,&Nrandom1);
  Node *random2Tree = readAndCreateTree(para,para.fileRanName2,para.ran2Id,rank,&Nrandom2);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(verbose) printf("D1D2...       ");  D1D2 = Npairs(&para,data1Tree,  data2Tree,   rank, size, firstCall,verbose);
  if(verbose) printf("D1R1...       ");  D1R1 = Npairs(&para,data1Tree,  random1Tree, rank, size, firstCall,verbose);
  if(verbose) printf("D2R2...       ");  D2R2 = Npairs(&para,data2Tree,  random2Tree, rank, size, firstCall,verbose);
  if(verbose) printf("R2R2...       ");  R1R2 = Npairs(&para,random1Tree,random2Tree, rank, size, firstCall,verbose);
  
  MPI_Send(D1D2,  para.nbins*(para.nboots+1), MPI_LONG, master, 1, MPI_COMM_WORLD);
  MPI_Send(D1R1,  para.nbins*(para.nboots+1), MPI_LONG, master, 2, MPI_COMM_WORLD);
  MPI_Send(D2R2,  para.nbins*(para.nboots+1), MPI_LONG, master, 3, MPI_COMM_WORLD);
  MPI_Send(R1R2,  para.nbins*(para.nboots+1), MPI_LONG, master, 4, MPI_COMM_WORLD);
  
  if(rank == master){
    long *D1D2Rank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    long *D1R1Rank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    long *D2R2Rank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    long *R1R2Rank = (long *)malloc(para.nbins*(para.nboots+1)*sizeof(long));
    
    for (k = 1; k < size; k++) {
      MPI_Recv(D1D2Rank, para.nbins*(para.nboots+1), MPI_LONG , k, 1, MPI_COMM_WORLD, &status);  
      MPI_Recv(D1R1Rank, para.nbins*(para.nboots+1), MPI_LONG , k, 2, MPI_COMM_WORLD, &status);  
      MPI_Recv(D2R2Rank, para.nbins*(para.nboots+1), MPI_LONG , k, 3, MPI_COMM_WORLD, &status);  
      MPI_Recv(R1R2Rank, para.nbins*(para.nboots+1), MPI_LONG , k, 4, MPI_COMM_WORLD, &status);  
      for(i = 0;i<para.nbins;i++){
	for(j=0;j<para.nboots+1;j++){
	  D1D2[j+(para.nboots+1)*i] += D1D2Rank[j+(para.nboots+1)*i];
	  D1R1[j+(para.nboots+1)*i] += D1R1Rank[j+(para.nboots+1)*i];
	  D2R2[j+(para.nboots+1)*i] += D2R2Rank[j+(para.nboots+1)*i];
	  R1R2[j+(para.nboots+1)*i] += R1R2Rank[j+(para.nboots+1)*i];
	}
      }
    }
    free(D1D2Rank);
    free(D1R1Rank);
    free(D2R2Rank);
    free(R1R2Rank);
  }
  
  if(rank == master)
    switch(para.estimator){
    case LS:
      writeLandySzalay(para,D1D2,D1R1,D2R2,R1R2,Ndata1,Nrandom1,Ndata2,Nrandom2);
      break;
    case NAT:
      writeNatural(para,D1D2,R1R2,Ndata1,Nrandom1,Ndata2,Nrandom2);
      break;
    case HAM:
      writeHamilton(para,D1D2,D1R1,D2R2,R1R2,Ndata1,Nrandom1,Ndata2,Nrandom2);
      break;
    }
  
  free(D1D2); 
  free(D1R1); 
  free(D2R2); 
  free(R1R2); 
  free_Node(data1Tree);
  free_Node(random1Tree);  
  free_Node(data2Tree);
  free_Node(random2Tree);  
  
  return SUCCESS;
}


void writeLandySzalay(Config para, long *D1D2, long *D1R1, long *D2R2, long *R1R2, long Ndata1, long Nrandom1, long Ndata2, long Nrandom2){
  /*Computes y(x) and errors from bootstrap. Then writes x y yerr in para.fileOutName.*/
  int i,j,n;
  double yi,yj,yn;
  
  //Initialization
  double Norm1 = (double)(Nrandom1*(Nrandom2-1))/(double)(Ndata1*(Ndata2-1));
  double Norm2 = (double)(Nrandom2-1)/(double)Ndata1;
  double Norm3 = (double)(Nrandom1-1)/(double)Ndata2;
  double *x    = (double *)malloc(para.nbins*sizeof(double));
  double *y    = (double *)malloc(para.nbins*sizeof(double));
  double *yerr = (double *)malloc(para.nbins*sizeof(double));

  //x
  for(i=0;i<para.nbins;i++){
    x[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
    if(para.log) x[i] = exp(x[i]);
  }
  
  //y
  for(i=0;i<para.nbins;i++){
    y[i]  =  Norm1*(double)D1D2[(para.nboots+1)*i]/(double)R1R2[(para.nboots+1)*i];
    y[i] += -Norm2*(double)D1R1[(para.nboots+1)*i]/(double)R1R2[(para.nboots+1)*i];
    y[i] += -Norm3*(double)D2R2[(para.nboots+1)*i]/(double)R1R2[(para.nboots+1)*i] + 1.0;
  }
  
  //Werr(theta)
  for(i=0;i<para.nbins;i++){
    yerr[i] = 0.0;
    for(n=1;n<para.nboots+1;n++){
      yn       =  Norm1*(double)D1D2[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i];
      yn      += -Norm2*(double)D1R1[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i];
      yn      += -Norm3*(double)D2R2[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i] + 1.0;
      yerr[i] += (y[i] - yn)*(y[i] - yn);
    }
    yerr[i] = sqrt(yerr[i]/(double)(para.nboots-1));
    //Computing half of the pairs (see Npairs()) overestimates 
    //bootstrap errors. A factor sqrt(2) is needed to correct for it 
    if(para.corr == AUTO) yerr[i] /= sqrt(2.0);
  }
  
  FILE *fileOut = fopen(para.fileOutName,"w");
  if(para.corr == AUTO){
    fprintf(fileOut,"#Auto-correlation. Landy & Szalay estimator\n");
    fprintf(fileOut,"#   coord     corr(coord)     corr_err           DD               DR            RR         Ndata      Nrandom\n");
    for(i=0;i<para.nbins;i++) fprintf(fileOut,"%12.7f %12.7f %12.7f %15zd %15zd %15zd %10zd %10zd\n",x[i],y[i],yerr[i],D1D2[(para.nboots+1)*i],D1R1[(para.nboots+1)*i],R1R2[(para.nboots+1)*i],Ndata1,Nrandom1);
  }else{
    fprintf(fileOut,"#Cross-correlation. Landy & Szalay estimator\n");
    fprintf(fileOut,"#   coord     corr(coord)     corr_err           D1D2             D1R1          D2R2         R1R2            Ndata1    Nrandom1   Ndata2    Nrandom2\n");
    for(i=0;i<para.nbins;i++) fprintf(fileOut,"%12.7f %12.7f %12.7f %15zd %15zd %15zd %15zd %10zd %10zd %10zd %10zd\n",x[i],y[i],yerr[i],D1D2[(para.nboots+1)*i],D1R1[(para.nboots+1)*i],D2R2[(para.nboots+1)*i],R1R2[(para.nboots+1)*i],Ndata1,Nrandom1,Ndata2,Nrandom2);
  }
  fclose(fileOut);
  
  if(para.cov_mat){
    FILE *fileCovOut = fopen(para.fileCovOutName,"w");
    double *ycov = (double *)malloc(para.nbins*para.nbins*sizeof(double));
    //Covariance matrix
    for(i=0;i<para.nbins;i++){
      for(j=0;j<para.nbins;j++){
	ycov[j+para.nbins*i] = 0.0;
	for(n=1;n<para.nboots+1;n++){
	  yi       =  Norm1*(double)D1D2[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i];
	  yi      += -Norm2*(double)D1R1[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i];
	  yi      += -Norm3*(double)D2R2[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i] + 1.0;
	  yj       =  Norm1*(double)D1D2[n+(para.nboots+1)*j]/(double)R1R2[n+(para.nboots+1)*j];
	  yj      += -Norm2*(double)D1R1[n+(para.nboots+1)*j]/(double)R1R2[n+(para.nboots+1)*j];
	  yj      += -Norm3*(double)D2R2[n+(para.nboots+1)*j]/(double)R1R2[n+(para.nboots+1)*j] + 1.0;
	  ycov[j+para.nbins*i] += (y[i] - yi)*(y[j] - yj)/(double)(para.nboots-1);
	}
	//Computing half of the pairs (see Npairs()) overestimates 
	//bootstrap errors. A factor 2 is needed to correct for it 
	if(para.corr == AUTO) ycov[j+para.nbins*i] /= 2.0;
	fprintf(fileCovOut,"%e   ", ycov[j+para.nbins*i]);
      }
      fprintf(fileCovOut,"\n");
    }
    fclose(fileCovOut);
  }
  
  return;
}

void writeNatural(Config para, long *D1D2, long *R1R2, long Ndata1, long Nrandom1, long Ndata2, long Nrandom2){
  /*Computes y(x) and errors from bootstrap. Then writes x y yerr in para.fileOutName.*/
  int i,j,n;
  double yi,yj,yn;
  
  //Initialization
  double Norm  = (double)(Nrandom1*(Nrandom2-1))/(double)(Ndata1*(Ndata2-1));
  double *x    = (double *)malloc(para.nbins*sizeof(double));
  double *y    = (double *)malloc(para.nbins*sizeof(double));
  double *yerr = (double *)malloc(para.nbins*sizeof(double));
  
  //x
  for(i=0;i<para.nbins;i++){
    x[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
    if(para.log) x[i] = exp(x[i]);
  }
  
  //y
  for(i=0;i<para.nbins;i++){
    y[i]  = Norm*(double)D1D2[(para.nboots+1)*i]/(double)R1R2[(para.nboots+1)*i] - 1.0;
  }
  
  //Werr(theta)
  for(i=0;i<para.nbins;i++){
    yerr[i] = 0.0;
    for(n=1;n<para.nboots+1;n++){
      yn       =  Norm*(double)D1D2[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i] - 1.0;
      yerr[i] += (y[i] - yn)*(y[i] - yn);
    }
    yerr[i] = sqrt(yerr[i]/(double)(para.nboots-1));
    //Computing half of the pairs (see Npairs()) overestimates 
    //bootstrap errors. A factor sqrt(2) is needed to correct for it 
    if(para.corr == AUTO) yerr[i] /= sqrt(2.0);
  }
  
  FILE *fileOut = fopen(para.fileOutName,"w");
  if(para.corr == AUTO){
    fprintf(fileOut,"#Auto-correlation. Natural estimator\n");
    fprintf(fileOut,"#   coord     corr(coord)     corr_err           DD                   RR         Ndata      Nrandom\n");
    for(i=0;i<para.nbins;i++) fprintf(fileOut,"%12.7f %12.7f %12.7f %15zd %15zd %10zd %10zd\n",x[i],y[i],yerr[i],D1D2[(para.nboots+1)*i],R1R2[(para.nboots+1)*i],Ndata1,Nrandom1);
  }else{
    fprintf(fileOut,"#Cross-correlation. Natural estimator\n");
    fprintf(fileOut,"#   coord     corr(coord)     corr_err           D1D2                 R1R2            Ndata1    Nrandom1   Ndata2    Nrandom2\n");
    for(i=0;i<para.nbins;i++) fprintf(fileOut,"%12.7f %12.7f %12.7f %15zd %15zd %10zd %10zd %10zd %10zd\n",x[i],y[i],yerr[i],D1D2[(para.nboots+1)*i],R1R2[(para.nboots+1)*i],Ndata1,Nrandom1,Ndata2,Nrandom2);
  }
  fclose(fileOut);
  
  if(para.cov_mat){
    FILE *fileCovOut = fopen(para.fileCovOutName,"w");
    double *ycov = (double *)malloc(para.nbins*para.nbins*sizeof(double));
    //Covariance matrix
    for(i=0;i<para.nbins;i++){
      for(j=0;j<para.nbins;j++){
	ycov[j+para.nbins*i] = 0.0;
	for(n=1;n<para.nboots+1;n++){
	  yi       =  Norm*(double)D1D2[n+(para.nboots+1)*i]/(double)R1R2[n+(para.nboots+1)*i] - 1.0;
	  yj       =  Norm*(double)D1D2[n+(para.nboots+1)*j]/(double)R1R2[n+(para.nboots+1)*j] - 1.0;
	  ycov[j+para.nbins*i] += (y[i] - yi)*(y[j] - yj)/(double)(para.nboots-1);
	}
	//Computing half of the pairs (see Npairs()) overestimates 
	//bootstrap errors. A factor 2 is needed to correct for it 
	if(para.corr == AUTO) ycov[j+para.nbins*i] /= 2.0;
	fprintf(fileCovOut,"%e   ", ycov[j+para.nbins*i]);
      }
      fprintf(fileCovOut,"\n");
    }
    fclose(fileCovOut);
  }
  return;
}


void writeHamilton(Config para, long *D1D2, long *D1R1, long *D2R2, long *R1R2, long Ndata1, long Nrandom1, long Ndata2, long Nrandom2){
  /*Computes y(x) and errors from bootstrap. Then writes x y yerr in para.fileOutName.*/
  int i,j,n;
  double yi,yj,yn;
  
  //Initialization
  //**Normalisation to be checked for cross correlation**
  double Norm  = (double)(Ndata2*Nrandom2)/(double)((Nrandom2-1)*(Ndata2-1));
  double *x    = (double *)malloc(para.nbins*sizeof(double));
  double *y    = (double *)malloc(para.nbins*sizeof(double));
  double *yerr = (double *)malloc(para.nbins*sizeof(double));

  //x
  for(i=0;i<para.nbins;i++){
    x[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
    if(para.log) x[i] = exp(x[i]);
  }
  
  //y
  for(i=0;i<para.nbins;i++){
    y[i] = Norm*(double)D1D2[(para.nboots+1)*i]*(double)R1R2[(para.nboots+1)*i]/((double)D1R1[(para.nboots+1)*i]*(double)D2R2[(para.nboots+1)*i]) - 1.0;
  }
  
  //Werr(theta)
  for(i=0;i<para.nbins;i++){
    yerr[i] = 0.0;
    for(n=1;n<para.nboots+1;n++){
      yn = Norm*(double)D1D2[n+(para.nboots+1)*i]*(double)R1R2[n+(para.nboots+1)*i]/((double)D1R1[n+(para.nboots+1)*i]*(double)D2R2[n+(para.nboots+1)*i]) - 1.0;
      yerr[i] += (y[i] - yn)*(y[i] - yn);
    }
    yerr[i] = sqrt(yerr[i]/(double)(para.nboots-1));
    //Computing half of the pairs (see Npairs()) overestimates 
    //bootstrap errors. A factor sqrt(2) is needed to correct for it 
    if(para.corr == AUTO) yerr[i] /= sqrt(2.0);
  }
  
  FILE *fileOut = fopen(para.fileOutName,"w");
  if(para.corr == AUTO){
    fprintf(fileOut,"#Auto-correlation. Hamilton estimator\n");
    fprintf(fileOut,"#   coord     corr(coord)     corr_err           DD               DR            RR         Ndata      Nrandom\n");
    for(i=0;i<para.nbins;i++) fprintf(fileOut,"%12.7f %12.7f %12.7f %15zd %15zd %15zd %10zd %10zd\n",x[i],y[i],yerr[i],D1D2[(para.nboots+1)*i],D1R1[(para.nboots+1)*i],R1R2[(para.nboots+1)*i],Ndata1,Nrandom1);
  }else{
    fprintf(fileOut,"#Cross-correlation. Hamilton estimator\n");
    fprintf(fileOut,"#   coord     corr(coord)     corr_err           D1D2             D1R1          D2R2         R1R2            Ndata1    Nrandom1   Ndata2    Nrandom2\n");
    for(i=0;i<para.nbins;i++) fprintf(fileOut,"%12.7f %12.7f %12.7f %15zd %15zd %15zd %15zd %10zd %10zd %10zd %10zd\n",x[i],y[i],yerr[i],D1D2[(para.nboots+1)*i],D1R1[(para.nboots+1)*i],D2R2[(para.nboots+1)*i],R1R2[(para.nboots+1)*i],Ndata1,Nrandom1,Ndata2,Nrandom2);
  }
  fclose(fileOut);
  
  if(para.cov_mat){
    FILE *fileCovOut = fopen(para.fileCovOutName,"w");
    double *ycov = (double *)malloc(para.nbins*para.nbins*sizeof(double));
    //Covariance matrix
    for(i=0;i<para.nbins;i++){
      for(j=0;j<para.nbins;j++){
	ycov[j+para.nbins*i] = 0.0;
	for(n=1;n<para.nboots+1;n++){
	  yi = Norm*(double)D1D2[n+(para.nboots+1)*i]*(double)R1R2[n+(para.nboots+1)*i]/((double)D1R1[n+(para.nboots+1)*i]*(double)D2R2[n+(para.nboots+1)*i]) - 1.0;
	  yj = Norm*(double)D1D2[n+(para.nboots+1)*j]*(double)R1R2[n+(para.nboots+1)*j]/((double)D1R1[n+(para.nboots+1)*j]*(double)D2R2[n+(para.nboots+1)*j]) - 1.0;
	  ycov[j+para.nbins*i] += (y[i] - yi)*(y[j] - yj)/(double)(para.nboots-1);
	}
	//Computing half of the pairs (see Npairs()) overestimates 
	//bootstrap errors. A factor 2 is needed to correct for it 
	if(para.corr == AUTO) ycov[j+para.nbins*i] /= 2.0;
	fprintf(fileCovOut,"%e   ", ycov[j+para.nbins*i]);
      }
      fprintf(fileCovOut,"\n");
    }
    fclose(fileCovOut);
  }
  
  return;
}

Node *readAndCreateTree(Config para,char *fileName, int ids[NIDS], int rank, long *N){
  int master = 0;
  FILE *fileIn;
  
  if(rank == master){
    fileIn = fopenAndCheck(fileName,"r");
    *N     = getCatSize(fileIn);
    //Sanity check ------------------------------
    if(*N < 1){
      printf("%s: input file has no point. Exiting...\n",MYNAME);
      exit(EXIT_FAILURE);
    }
  }
  MPI_Bcast (N,1,MPI_LONG,master,MPI_COMM_WORLD);
  double *data = (double *)malloc(NIDS*(*N)*sizeof(double));
  if(rank == master){
    printf("Reading %s...",fileName);fflush(stdout);
    readCat(fileIn,ids,data);
    printf("(%zd objects found). Building tree...",*N);fflush(stdout);
    fclose(fileIn);
  }
  MPI_Bcast (data,NIDS*(*N),MPI_DOUBLE,master,MPI_COMM_WORLD);
  Point *dataPoint = dataToPoint(para,data,*N);
  Node *dataTree   = createNode(para,dataPoint,*N,0,1);
  
  if(rank == master) printf("\n");
  
  free(data);
  free_Point(dataPoint,*N);
  return dataTree;
} 

/*----------------------------------------------------------------*
 *Tree routines                                                   *
 *----------------------------------------------------------------*/

Node  **getNodesForCpu(Config *para, Node *data, int count, int firstCall){
  /*Splits the tree into subtrees to be treated separately by individual 
    cpus and returns the array of nodes. "count" should be the 
    total number of cpus.*/
  
  static int index;
  static Node **result;
  
  if(firstCall){
    index = 0;
    result = (Node **)malloc(count*sizeof(Node));
  }
  
  if(count > 2 && count%2 > 0){
    getNodesForCpu(para,data->Left,(count-1)/2,0);
    getNodesForCpu(para,data->Right,(count+1)/2,0);
  }else if(count > 2 && count%2==0){
    getNodesForCpu(para,data->Left,count/2,0);
    getNodesForCpu(para,data->Right,count/2,0);
  }else if(count == 2){
    result[index] = (Node *)data->Left;
    index++;
    result[index] = (Node *)data->Right;
    index++;
    return result;
  }else{
    result[index] = (Node *)data;
    index++;
    return result;
  }  
  return result;
}


void gg(Config *para, Node *lens, Node *source, double *SigR, double *weight, double *info, int rank, int size,  int firstCall, int verbose){
  /* Computes the galaxy-galaxy two-point correlation function */
  double d;
  Node **node1Rank;
  static long i,total, count;  
  
  if(firstCall){
    /* distribute tree nodes among cpus */
    node1Rank = getNodesForCpu(para,source,size,firstCall);
    source    = node1Rank[rank];
    for(i=0;i<size;i++){
      if(i != rank) free_Node(node1Rank[i]);
    }
    total = lens->N[0]*source->N[0]; 
    count = 0;
    fflush(stdout);
  }
  d = para->distAng(&lens->point[0],&source->point[0]);
  
  /* correlate nodes if (type = LEAF OR size/d < OA) AND Deltaz < sigmaz*(1+z)) */
  if(lens->type != LEAF && lens->radius/d > para->OA || lens->Deltaz > para->sigz*(1.0+lens->point[0].x[2])){
    if( source->type != LEAF && source->radius/d > para->OA || source->Deltaz > para->sigz*(1.0+source->point[0].x[2])){
      gg(para,lens->Left,  source->Left,  SigR, weight, info, rank, size, 0, verbose);
      gg(para,lens->Right, source->Right, SigR, weight, info, rank, size, 0, verbose);
      gg(para,lens->Left,  source->Right, SigR, weight, info, rank, size, 0, verbose);
      gg(para,lens->Right, source->Left,  SigR, weight, info, rank, size, 0, verbose);
    }else{
      gg(para,lens->Left,source,  SigR, weight, info, rank, size, 0, verbose);
      gg(para,lens->Right,source, SigR, weight, info, rank, size, 0, verbose);
    }
  }else if( source->type != LEAF && source->radius/d > para->OA || source->Deltaz > para->sigz*(1.0+source->point[0].x[2])){
    gg(para,lens,source->Left,  SigR, weight, info, rank, size, 0, verbose);
    gg(para,lens,source->Right, SigR, weight, info, rank, size, 0, verbose);
  }else{
    corrLensSource(para,lens->point[0], source->point[0], d, SigR, weight, info);
    count +=  lens->N[0]*source->N[0];
    printCount(&count, &total, verbose);
  }
  
  if(firstCall){
    free_Node(source);
    free_Node(lens);
    if(verbose) printf("\b\b\b\b\b\b\b%6.2f%%\n",100.0);
  }
  return;
}

void corrLensSource(Config *para, Point lens, Point source, double d, double *SigR, double *weight, double *info){
  /* See Leauthaud et al. (2010),  2010ApJ...709...97L */
  double dA, dR, AS, AL;
  double phi_gg, cos2phi_gg, sin2phi_gg, e1, e2, w;
  double fac, DOS, DOL, DLS, SigCritInv;
  Point *A = createPoint(1,para->nboots+1);
  
  int i,n,k;
  
  switch(para->proj){
  case COMO:
    fac = lens.x[2]+1.0; /* (1+z_lens)*/
    break;
  case PHYS:
    fac = 1.0;
    break;
  case THETA:
    printf("Wrong coordinate system for gal-gal lensing\n");
    exit(EXIT_FAILURE);
    break;
  }
  
  /* DA = phy_dis/angular_size_in_radians = tran_como_dis/(1+z)*/
  
  dA  = dComo(lens.x[2],para->a)/(lens.x[2]+1.0); /* Angular diameter distance in comoving coordinates*/ 
  dR  = d*dA*PI/180.0;                            /* Transverse distance in phys coordinates (Mpc)    */
  dR *= fac;                                      /* If coordinates in comoving system                */
  if(para->log) dR = log(dR);
  k = floor((dR-para->min)/para->Delta);
  if(0 <= k && k < para->nbins && source.x[2] > lens.x[2] + source.sigz + lens.sigz){
    //Point A ---------------------------------
    A->x[0]    = source.x[0];
    A->cosx[0] = source.cosx[0];
    A->sinx[0] = source.sinx[0];
    A->x[1]    = lens.x[1];
    A->cosx[1] = lens.cosx[1];
    A->sinx[1] = lens.sinx[1];
    AS         = para->distAng(A,&source);
    AL         = para->distAng(A,&lens);
    //To get correct sign for phi_gg -----------
    if(source.x[0] > lens.x[0]) AL = -AL;
    if(source.x[1] < lens.x[1]) AS = -AS;
    phi_gg     = atan2(AS,AL);
    cos2phi_gg = cos(2.0*phi_gg);
    sin2phi_gg = sin(2.0*phi_gg);
    e1         =  source.e1*cos2phi_gg + source.e2*sin2phi_gg;
    e2         = -source.e1*sin2phi_gg + source.e2*cos2phi_gg;
    //Lens equation -----------------------------
    DOS        = dComo(source.x[2],para->a);/* Note this is a comoving los distance,         */ 
    DOL        = dComo(lens.x[2],para->a);  /* but redshift factors cancel out               */
    DLS        = (dComo(source.x[2],para->a)-dComo(lens.x[2],para->a));
    SigCritInv = DOL*DLS/DOS/1000.0/fac;           /* 1/SigCrit in Gpc                       */
    SigCritInv /= 1.663e3;                         /* see Narayan & Bartelman pge 10         */
    SigCritInv *= fac*fac;                         /* If coordinates in comoving system      */
    w = SigCritInv*SigCritInv*source.w;
    for(n=0;n<para->nboots+1;n++){
      weight[n+(para->nboots+1)*k]  += lens.boots[n]*w;
      SigR[n+(para->nboots+1)*k]    += lens.boots[n]*e1*w/SigCritInv;
    }
    /* keep track of info per bin in this memory slot: */
    info[0+para->Ninfo*k] += 1.0;             /* Nsources */
    info[1+para->Ninfo*k] += dR*w;            /* mean R */
    info[2+para->Ninfo*k] += e2*w/SigCritInv; /* e2       */
  }
  free_Point(A,1);
}

long *Npairs(Config *para, Node *node1, Node *node2, int rank, int size,  int firstCall, int verbose){
  double d, dA;
  int i,j,k,l;
  static long *NN, total, count, fac;
  
  if(firstCall){
    /* split tree among cpus */
    Node **node1Rank   = getNodesForCpu(para,node1,size,firstCall);
    node1 = node1Rank[rank];
    /* initialization */
    NN = (long *)malloc(para->nbins*(para->nboots+1)*sizeof(long));
    for(i=0;i<para->nbins;i++){
       for(j=0;j<para->nboots+1;j++) NN[j+(para->nboots+1)*i] = 0;
    }
    total = node1->N[0]*node2->N[0]; 
    count = 0;
    fflush(stdout);
    /* if node1 and node2 are the same, only compute half of pairs */
    fac = 1;
    if(node1->root == node2->root) fac = 2;
  }
  
  if(node1->root == node2->root && node1->id < node2->id){
    count +=  node1->N[0]*node2->N[0];
    printCount(&count,&total,verbose);
    return NN;
  }
  
  d = para->distAng(&node1->point[0],&node2->point[0]);
  
  if(node1->type != LEAF && node1->radius/d > para->OA){
    if( node2->type != LEAF && node2->radius/d > para->OA){
      Npairs(para,node1->Left,  node2->Left,  rank, size, 0, verbose);
      Npairs(para,node1->Right, node2->Right, rank, size, 0, verbose);
      Npairs(para,node1->Left,  node2->Right, rank, size, 0, verbose);
      Npairs(para,node1->Right, node2->Left,  rank, size, 0, verbose);
    }else{
      Npairs(para,node1->Left,node2,  rank, size, 0, verbose);
      Npairs(para,node1->Right,node2, rank, size, 0, verbose);
    }
  }else if( node2->type != LEAF && node2->radius/d > para->OA){
    Npairs(para,node1,node2->Left,  rank, size, 0, verbose);
    Npairs(para,node1,node2->Right, rank, size, 0, verbose);
  }else{
    if(para->proj == PHYS){
      /* angular diameter distance in comoving coordinates */ 
      dA = dComo(node1->point[0].x[2],para->a)/(node1->point[0].x[2]+1.0); 
      /* transverse distance in phys coordinates (Mpc)     */
      d  = d*dA*PI/180.0;                           
    }
    if(para->log) d = log(d);
    k = floor((d-para->min)/para->Delta);
    if(0 <= k && k < para->nbins){
      for(i=0;i<para->nboots+1;i++) NN[i+(para->nboots+1)*k] += fac*node1->N[i]*node2->N[i];
    }
    count +=  node1->N[0]*node2->N[0];
    printCount(&count,&total,verbose);
  }
  if(firstCall && verbose) printf("\b\b\b\b\b\b\b%6.2f%%\n",100.0);
  
  return NN;
}



Node *createNode2(Config para, Point *data, long N, int SplitDim, int firstCall){

  //THIS IS TO TEST METHODS AIMED AT REDUCING THE MEMORY USAGE
  


  /*Returns a pointer to the node, recursively. 
    If N < Nmin it returns a node with no child and sets node.Left = NULL 
    and node.Right = NULL (:= leaf). "N" is the number of points INSIDE the node.*/
  long i,j,Nmin = 1;
  static void *root;
  static long countNodes;
  static Point *dataRoot, *workspace, *save;
  static double *xmin,*xmax;
  
  Node *result  = (Node *)malloc(sizeof(Node));    /* Allocate memory for THIS node */
  
  if(firstCall){                                   /* Allocate workspaces */ 
    xmin       = (double *)malloc(NDIM*sizeof(double));
    xmax       = (double *)malloc(NDIM*sizeof(double));
    root       = result;
    countNodes = 0;
    workspace  = createPoint(N,para.nboots+1);
    save       = createPoint(N,para.nboots+1); 
    for(i=0;i<N;i++) cpyPoints(&save[i],&data[i]); /* Save data */
  }
  result->root = root;
  result->id   = countNodes;
  countNodes++;  
  
  if(N <= Nmin){               /* leaf */
    result->Npoint = 1;
    result->point  =  createPoint(result->Npoint,para.nboots+1);
    result->N      = (long *)malloc((para.nboots+1)*sizeof(long));
    for(i=0;i<para.nboots+1;i++) result->N[i] = data[0].boots[i];
    cpyPoints(&result->point[0],&data[0]);
    result->radius = 0.0;
    result->Deltaz = 0.0;
    result->type   = LEAF;
    result->Left   = NULL;
    result->Right  = NULL;
  }else{                        /* node */
    result->type = NODE;
    for(j=0;j<NDIM;j++){        /* node limits */
      xmin[j] = xmax[j] = data[0].x[j];
      for(i=0;i<N;i++) {
	xmin[j] = MIN(xmin[j],data[i].x[j]);
	xmax[j] = MAX(xmax[j],data[i].x[j]);
      }
    }
    if(NDIM > 2) result->Deltaz = xmax[2] - xmin[2];

    
   
    result->Npoint = 1;                       
    result->point  =  createPoint(result->Npoint,para.nboots+1);
    
    result->point[0].x[0] = (xmax[0] + xmin[0])/2.0;
    result->point[0].x[1] = (xmax[1] + xmin[1])/2.0;
    result->point[0].cosx[0] = cos(result->point[0].x[0]*PI/180.0);
    result->point[0].sinx[0] = sin(result->point[0].x[0]*PI/180.0);
    result->point[0].cosx[1] = cos(result->point[0].x[1]*PI/180.0);
    result->point[0].sinx[1] = sin(result->point[0].x[1]*PI/180.0);

  
      //Node size
    result->radius = EPS;
    for(i=0;i<N;i++){
      result->radius = MAX(result->radius,para.distAng(&result->point[0], &data[i]));
    }
    
    //  if(result->radius > para.OA*((Node *)root)->radius){

    //  result->N      = (long *)malloc(1*sizeof(long));  /* number of objects */
    //  result->N[0]   = N;

    //    }else{
      result->N      = (long *)malloc((para.nboots+1)*sizeof(long));  /* number of objects */
      for(i=0;i<para.nboots+1;i++){
	result->N[i] = 0;
	for(j=0;j<N;j++) result->N[i] += data[j].boots[i];
      }
      for(j=0;j<NDIM;j++){
	result->point[0].x[j] = 0.0;
	for(i=0;i<N;i++) {
	  result->point[0].x[j] += data[i].x[j]/(double)N;
	}
	if(j<2){
	  result->point[0].cosx[j] = cos(result->point[0].x[j]*PI/180.0);
	  result->point[0].sinx[j] = sin(result->point[0].x[j]*PI/180.0);
	}
      }  
      result->point[0].boots[0] = result->N[0];
      for(i=1;i<para.nboots+1;i++) result->point[0].boots[i] = roundToNi((double)result->N[i]/(double)N); 
      if(para.corr == GGLENS){
	result->point[0].sigz = 0.0;
	result->point[0].e1   = 0.0;
	result->point[0].e2   = 0.0;
	result->point[0].w    = 0.0;
	for(i=0;i<N;i++) {
	  result->point[0].sigz += data[i].sigz/(double)N;
	  result->point[0].e1   += data[i].e1/(double)N;
	  result->point[0].e2   += data[i].e2/(double)N;
	  result->point[0].w    += data[i].w/(double)N;
	}
      }
      //  }
    
    
    double SplitValue;
    long NLeft, NRight;
    
    for(i=0;i<N;i++){   /* sort values along SplitDim coordinates */
      cpyPoints(&workspace[i], &data[i]);
      workspace[i].SplitDim = SplitDim;
    }
    qsort(workspace,N,sizeof(Point),comparePoints);
    
    //Find SplitValue
    switch(PARITY(N)){
    case EVEN:
      NLeft  = N/2;
      NRight = N/2;
      break;
    case ODD:
      NLeft  = (N+1)/2;
      NRight = (N-1)/2;
      break;
    }
    
    //Next splitting coordinate
    SplitDim++;
    if(SplitDim > NDIM-1)  SplitDim = 0;
  
    //Left ----------------------------------------------------
    for(i=0;i<NLeft;i++) cpyPoints(&data[i], &workspace[i]);
    result->Left = createNode(para,data,NLeft,SplitDim,0);
    
    //Right ----------------------------------------------------
   for(i=NLeft;i<N;i++) cpyPoints(&data[i - NLeft], &workspace[i]);
    result->Right = createNode(para,data,NRight,SplitDim,0);
  }
  
  if(firstCall){
    free(xmin);
    free(xmax);
    free_Point(workspace,N);
    for(i=0;i<N;i++) cpyPoints(&data[i],&save[i]); /* Restore data */
  }
  return result;
}

Node *createNode(Config para, Point *data, long N, int SplitDim, int firstCall){
  /*Returns a pointer to the node, recursively. 
    If N < Nmin it returns a node with no child and sets node.Left = NULL 
    and node.Right = NULL (:= leaf). "N" is the number of points INSIDE the node.*/
  long i,j,Nmin = 1;
  static void *root;
  static long countNodes;
  static Point *dataRoot, *workspace, *save;
  double zmin,zmax;
  
  Node *result  = (Node *)malloc(sizeof(Node)); /* Allocate memory for THIS node */
  
  if(firstCall){                     
    root       = result;
    countNodes = 0;
    workspace  = createPoint(N,para.nboots+1);
    save       = createPoint(N,para.nboots+1); 
    for(i=0;i<N;i++) cpyPoints(&save[i],&data[i]); /* Save data */
  }
  result->root = root;
  result->id   = countNodes;
  countNodes++;
  
  //Number of objects
  result->N      = (long *)malloc((para.nboots+1)*sizeof(long));
  for(i=0;i<para.nboots+1;i++){
    result->N[i] = 0;
    for(j=0;j<N;j++) result->N[i] += data[j].boots[i];
  }
  
  //Point with mean properties of the node
  result->Npoint = 1;
  result->point  =  createPoint(result->Npoint,para.nboots+1);
  for(j=0;j<NDIM;j++){
    result->point[0].x[j] = 0.0;
    for(i=0;i<N;i++) {
      result->point[0].x[j] += data[i].x[j]/(double)N;
    }
    if(j<2){
      result->point[0].cosx[j] = cos(result->point[0].x[j]*PI/180.0);
      result->point[0].sinx[j] = sin(result->point[0].x[j]*PI/180.0);
    }
  }  
  result->point[0].boots[0] = result->N[0];
  for(i=1;i<para.nboots+1;i++) result->point[0].boots[i] = roundToNi((double)result->N[i]/(double)N); 
  if(para.corr == GGLENS){
    result->point[0].sigz = 0.0;
    result->point[0].e1   = 0.0;
    result->point[0].e2   = 0.0;
    result->point[0].w    = 0.0;
    zmin = 99.0;
    zmax = -99.0;
    for(i=0;i<N;i++) {
      result->point[0].sigz += data[i].sigz/(double)N;
      result->point[0].e1   += data[i].e1/(double)N;
      result->point[0].e2   += data[i].e2/(double)N;
      result->point[0].w    += data[i].w/(double)N;
      zmin = MIN(zmin,data[i].x[2]);
      zmax = MAX(zmax,data[i].x[2]);
    }
    result->Deltaz = zmax - zmin;
  }
  
  //Node size
  result->radius = EPS;
  for(i=0;i<N;i++){
    result->radius = MAX(result->radius,para.distAng(&result->point[0], &data[i]));
  }
  
  if(N <= Nmin){
    //Leaf -----------------------------------------------------------------//
    result->type   = LEAF;
    result->Left   = NULL;
    result->Right  = NULL;
  }else{
    //Node -----------------------------------------------------------------//
    result->type = NODE;
    
    double SplitValue;
    long NLeft, NRight;
    
    for(i=0;i<N;i++){   /* sort values along SplitDim coordinates */
      cpyPoints(&workspace[i], &data[i]);
      workspace[i].SplitDim = SplitDim;
    }
    qsort(workspace,N,sizeof(Point),comparePoints);
    
    //Find SplitValue
    switch(PARITY(N)){
    case EVEN:
      NLeft  = N/2;
      NRight = N/2;
      break;
    case ODD:
      NLeft  = (N+1)/2;
      NRight = (N-1)/2;
      break;
    }
    
    //Next splitting coordinate
    SplitDim++;
    if(SplitDim > NDIM-1)  SplitDim = 0;
  
    //Left ----------------------------------------------------
    for(i=0;i<NLeft;i++) cpyPoints(&data[i], &workspace[i]);
    result->Left = createNode(para,data,NLeft,SplitDim,0);
    
    //Right ----------------------------------------------------
   for(i=NLeft;i<N;i++) cpyPoints(&data[i - NLeft], &workspace[i]);
    result->Right = createNode(para,data,NRight,SplitDim,0);
  }
  
  if(firstCall){
    free_Point(workspace,N);
    for(i=0;i<N;i++) cpyPoints(&data[i],&save[i]); /* Restore data */
  }
  return result;
}

/*----------------------------------------------------------------*
 *Configuration                                                   *
 *----------------------------------------------------------------*/

int readPara(int argc, char **argv, Config *para){
  /*Reads the input parameters given in the configuration file and returns
    para, that contains all the configuration parameters.*/
  int i,j,noconfigFile=1;
  long Ncol;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR], list[NFIELD*NCHAR], *str_end;
  
  //Initialisation --------------------------------
  srand((unsigned int)time(NULL));
  EPS   = determineMachineEpsilon();
  IDERR = determineIntError();
  strcpy(MYNAME,argv[0]);
  
  //Default parameters ----------------------------
  int ids[NIDS] = {1,2,3,4,5,6,7};
  strcpy(para->fileInName1,"data1.cat");
  strcpy(para->fileInName2,"data2.cat");
  strcpy(para->fileRanName1,"ran1.cat");
  strcpy(para->fileRanName2,"ran2.cat");
  for(j=0;j<NIDS;j++){
    para->data1Id[j] = ids[j];
    para->data2Id[j] = ids[j];
    para->ran1Id[j]  = ids[j];
    para->ran2Id[j]  = ids[j];
  }
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
  para->nboots    = 50;
  para->cov_mat   = 0;
  strcpy(para->fileOutName,   "corr.out");
  strcpy(para->fileCovOutName,"cov.out");
  
  //Cosmology
  para->a[0] = H0;
  para->a[1] = Omega_M;
  para->a[2] = Omega_L;
  para->a[3] = c;
  para->sigz = 0.03;
  
  for(i=0;i<argc;i++){
    //Help-------------------------------------------------------------------//
    if(!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help") || argc == 1){
      printf("\n\n\
                          S W O T\n\n\
                (Super W Of Theta) MPI version 1.42\n\n\
              (LEGACY VERSION - DO NOT USE)\n\n\
              (LEGACY VERSION - DO NOT USE)\n\n\
              (LEGACY VERSION - DO NOT USE)\n\n\
Program to compute two-point correlation functions.\n\
Usage:  %s -c configFile [options]: run the program\n\
        %s -d: display a default configuration file\n\
Important: if using \"RADEC\" coordinates, the angle\n\
in the input catalogues must be in decimal degrees.\n\
Output format: coord corr(coord) corr_err\n",MYNAME,MYNAME);
      exit(EXIT_FAILURE);
    }
    //Create a default configuration file ----------------------------
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
      printf("coord          RADEC\t  #Coordinates: RADEC (in deg) or CART\n");
      printf("#----------------------------------------------------------#\n");
      printf("#Correlation options                                       #\n");
      printf("#----------------------------------------------------------#\n");
      printf("corr           auto\t #Type of correlation: auto, cross or gglens\n");
      printf("est            ls\t #Estimator: ls, nat, ham\n");
      printf("range          %g,%g\t #Correlation range\n",para->min,para->max);
      printf("nbins          %d\t #Number of bins\n",para->nbins);
      printf("log            yes\t #Logarithmic bins (yes or no)\n");
      printf("nboots         %d\t #Number of bootstrap samples\n",para->nboots);
      printf("OA             %g\t #Open angle for approximation\n",para->OA);
      printf("proj           theta\t #Axis projection (or phys)\n");
      printf("#----------------------------------------------------------#\n");
      printf("#Cosmology (for gal-gal correlations) and w(R)             #\n");
      printf("#----------------------------------------------------------#\n");
      printf("H0             %g\t #Hubble parameter\n",para->a[0]);
      printf("Omega_M        %g\t #Relative matter density\n",para->a[1]);
      printf("Omega_L        %g\t #Relative energy density (Lambda)\n",para->a[2]);
      printf("sigz           %g\t #Redshift accuracy\n",para->sigz);
      printf("#proj           como\t #Axis projection (or phys)\n");
      printf("#----------------------------------------------------------#\n");
      printf("#Output options                                            #\n");
      printf("#----------------------------------------------------------#\n");
      printf("out            %s\t #Output file: coord corr(coord) corr_err\n",para->fileOutName);
      printf("cov            no\t #Covariance matrix of bootstrap errors (yes or no)\n");
      printf("cov_out        %s\t #Covariance output file\n",para->fileCovOutName);
      exit(EXIT_FAILURE);
    }
    //Read configuration file ----------------------------------
    if(!strcmp(argv[i],"-c")){
      FILE *configFile = fopenAndCheck(argv[i+1],"r");
      noconfigFile     = 0;
      while(fgets(line,NFIELD*NCHAR,configFile) != NULL){
	if(getStrings(line,item," ",&Ncol)){
	  if(!strcmp(getCharValue(item,1),"data1")) strcpy(para->fileInName1,getCharValue(item,2));
	  if(!strcmp(getCharValue(item,1),"cols1")){ 
	    getStrings(getCharValue(item,2),list,",",&Ncol);
	    for(j=0;j<Ncol;j++)  para->data1Id[j] = getIntValue(list,j+1);
	  }
	  if(!strcmp(getCharValue(item,1),"data2")) strcpy(para->fileInName2,getCharValue(item,2));
	  if(!strcmp(getCharValue(item,1),"cols2")){ 
	    getStrings(getCharValue(item,2),list,",",&Ncol);
	    for(j=0;j<Ncol;j++)  para->data2Id[j] = getIntValue(list,j+1);
	  }
	  if(!strcmp(getCharValue(item,1),"ran1")) strcpy(para->fileRanName1,getCharValue(item,2));
	  if(!strcmp(getCharValue(item,1),"rancols1")){ 
	    getStrings(getCharValue(item,2),list,",",&Ncol);
	    for(j=0;j<Ncol;j++)  para->ran1Id[j] = getIntValue(list,j+1);
	  }
	  if(!strcmp(getCharValue(item,1),"ran2")) strcpy(para->fileRanName2,getCharValue(item,2));
	  if(!strcmp(getCharValue(item,1),"rancols2")){ 
	    getStrings(getCharValue(item,2),list,",",&Ncol);
	    for(j=0;j<Ncol;j++)  para->ran2Id[j] = getIntValue(list,j+1);
	  }
	  if(!strcmp(getCharValue(item,1),"coord")) {
	    if(!strcmp(getCharValue(item,2),"RADEC")) para->coordType = RADEC;
	    if(!strcmp(getCharValue(item,2),"CART"))  para->coordType = CART;
	  }
	  if(!strcmp(getCharValue(item,1),"corr")) {
	    if(!strcmp(getCharValue(item,2),"cross"))  para->corr = CROSS;
	    if(!strcmp(getCharValue(item,2),"auto"))   para->corr = AUTO;
	    if(!strcmp(getCharValue(item,2),"gglens")) para->corr = GGLENS;
	  }
	  if(!strcmp(getCharValue(item,1),"proj")) {
	    if(!strcmp(getCharValue(item,2),"theta"))  para->proj = THETA;
	    if(!strcmp(getCharValue(item,2),"como"))   para->proj = COMO;
	    if(!strcmp(getCharValue(item,2),"phys"))   para->proj = PHYS;
	  }
	  if(!strcmp(getCharValue(item,1),"est")) {
	    if(!strcmp(getCharValue(item,2),"ls"))  para->estimator = LS;
	    if(!strcmp(getCharValue(item,2),"nat")) para->estimator = NAT;
	    if(!strcmp(getCharValue(item,2),"ham")) para->estimator = HAM;
	  }
	  if(!strcmp(getCharValue(item,1),"range")){ 
	    getStrings(getCharValue(item,2),list,",",&Ncol);
	    para->min = getDoubleValue(list,1);
	    para->max = getDoubleValue(list,2);
	  }
	  if(!strcmp(getCharValue(item,1),"nbins")) para->nbins = getIntValue(item,2);
	  if(!strcmp(getCharValue(item,1),"log")){
	    if(!strcmp(getCharValue(item,2),"yes")) para->log = 1;
	    else para->log = 0;
	  }
	  if(!strcmp(getCharValue(item,1),"nboots"))  para->nboots = getIntValue(item,2);
	  if(!strcmp(getCharValue(item,1),"OA"))      para->OA     = getDoubleValue(item,2);
	  if(!strcmp(getCharValue(item,1),"H0"))      para->a[0]   = getDoubleValue(item,2);
	  if(!strcmp(getCharValue(item,1),"Omega_M")) para->a[1]   = getDoubleValue(item,2);
	  if(!strcmp(getCharValue(item,1),"Omega_L")) para->a[2]   = getDoubleValue(item,2);
	  if(!strcmp(getCharValue(item,1),"sigz"))    para->sigz   = getDoubleValue(item,2);
	  if(!strcmp(getCharValue(item,1),"out") || !strcmp(getCharValue(item,1),"o")) strcpy(para->fileOutName,getCharValue(item,2));
	  if(!strcmp(getCharValue(item,1),"cov")){
	    if(!strcmp(getCharValue(item,2),"yes")) para->cov_mat = 1;
	    else para->cov_mat = 0;
	  }
	  if(!strcmp(getCharValue(item,1),"cov_out")) strcpy(para->fileCovOutName,getCharValue(item,2));
	}
      }
      fclose(configFile);
    }
    //Read command line options ----------------------------------
    //getCharValue(item,1) -> argv[i],getCharValue(item,2) -> argv[i+1]
    if(!strcmp(argv[i],"-data1")) strcpy(para->fileInName1,argv[i+1]);
    if(!strcmp(argv[i],"-cols1")){ 
      getStrings(argv[i+1],list,",",&Ncol);
      for(j=0;j<Ncol;j++)  para->data1Id[j] = getIntValue(list,j+1);
    }
    if(!strcmp(argv[i],"-data2")) strcpy(para->fileInName2,argv[i+1]);
    if(!strcmp(argv[i],"-cols2")){ 
      getStrings(argv[i+1],list,",",&Ncol);
      for(j=0;j<Ncol;j++)  para->data2Id[j] = getIntValue(list,j+1);
    }
    if(!strcmp(argv[i],"-ran1")) strcpy(para->fileRanName1,argv[i+1]);
    if(!strcmp(argv[i],"-rancols1")){ 
      getStrings(argv[i+1],list,",",&Ncol);
      for(j=0;j<Ncol;j++)  para->ran1Id[j] = getIntValue(list,j+1);
    }
    if(!strcmp(argv[i],"-ran2")) strcpy(para->fileRanName2,argv[i+1]);
    if(!strcmp(argv[i],"-rancols2")){ 
      getStrings(argv[i+1],list,",",&Ncol);
      for(j=0;j<Ncol;j++)  para->ran2Id[j] = getIntValue(list,j+1);
    }
    if(!strcmp(argv[i],"-coord")) {
      if(!strcmp(argv[i+1],"RADEC")) para->coordType = RADEC;
      if(!strcmp(argv[i+1],"CART"))  para->coordType = CART;
    }
    if(!strcmp(argv[i],"-corr")) {
      if(!strcmp(argv[i+1],"cross"))  para->corr = CROSS;
      if(!strcmp(argv[i+1],"auto"))   para->corr = AUTO;
      if(!strcmp(argv[i+1],"gglens")) para->corr = GGLENS;
    }
    if(!strcmp(argv[i],"-proj")) {
      if(!strcmp(argv[i+1],"theta")) para->proj = THETA;
      if(!strcmp(argv[i+1],"como"))  para->proj = COMO;
      if(!strcmp(argv[i+1],"phys"))  para->proj = PHYS;
    }
    if(!strcmp(argv[i],"-est")) {
      if(!strcmp(argv[i+1],"ls"))  para->estimator = LS;
      if(!strcmp(argv[i+1],"nat")) para->estimator = NAT;
      if(!strcmp(argv[i+1],"ham")) para->estimator = HAM;
    }
    if(!strcmp(argv[i],"-range")){ 
      getStrings(argv[i+1],list,",",&Ncol);
      para->min = getDoubleValue(list,1);
      para->max = getDoubleValue(list,2);
    }
    if(!strcmp(argv[i],"-nbins")) para->nbins = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-log")){
      if(!strcmp(argv[i+1],"yes")) para->log = 1;
      else para->log = 0;
    }
    if(!strcmp(argv[i],"-nboots"))  para->nboots = atoi(argv[i+1]);
    if(!strcmp(argv[i],"-OA"))      para->OA     = atof(argv[i+1]);
    if(!strcmp(argv[i],"-H0"))      para->a[0]   = atof(argv[i+1]);
    if(!strcmp(argv[i],"-Omega_M")) para->a[1]   = atof(argv[i+1]);
    if(!strcmp(argv[i],"-Omega_L")) para->a[2]   = atof(argv[i+1]);
    if(!strcmp(argv[i],"-sigz"))    para->sigz   = atof(argv[i+1]);
    if(!strcmp(argv[i],"-out") || !strcmp(argv[i],"-o")) strcpy(para->fileOutName,argv[i+1]);
    if(!strcmp(argv[i],"-cov")){
      if(!strcmp(argv[i+1],"yes")) para->cov_mat = 1;
      else para->cov_mat = 0;
    }
    if(!strcmp(argv[i],"-cov_out")) strcpy(para->fileCovOutName,argv[i+1]);
  }
  
  if(noconfigFile){
    printf("Please provide a configuration file: swot -c configFile\n");
    exit(EXIT_FAILURE);
  }
  
  if(para->corr == GGLENS) NDIM  = 3;
  if(para->corr == GGLENS && para->proj == THETA)  para->proj = COMO; /*Default projection for gg lensing*/

  switch(para->coordType)
    {
    case RADEC:
      para->distAng = &distAngSpher;
      break;
    case CART:
      para->distAng = &distAngCart;
      break;
    }
  
  if(para->log){
    para->min = log(para->min);
    para->max = log(para->max);
  }  
  para->Delta = (para->max - para->min)/(double)para->nbins;
  
  return SUCCESS;
}

/*----------------------------------------------------------------*
 *Utils - Numeric                                                 *
 *----------------------------------------------------------------*/


double dComo(double z, double a[4]){
  /*Return the comoving distance of z*/

  int n = 1000;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (n);
  
  gsl_function F;
  F.function = &drdz;
  F.params   = a;

  double result,error;
  gsl_integration_qag(&F, 0.0, z, 0, 1e-7, n, 6, w, &result, &error); 
  gsl_integration_workspace_free (w);
  
  return result;
}

double drdz(double x, void * params){
  /*
    a[0] = H0;
    a[1] = Omega_M;
    a[2] = Omega_L;
    a[3] = c;
  */
  
  double *a = (double *) params;
  return a[3]/(a[0]*sqrt(a[1]*pow(1+x,3.0)+a[2]));
}

int getCatSize(FILE *fileIn){
  /*Returns the number of lines in fileIn*/
  long i,Ncol;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR];
  
  int result = 0;
  rewind(fileIn);
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL)
    if(*(line) != '#' && *(line) != '\0' && *(line) != '\n') result++;
  rewind(fileIn);
  return result;
}

void readCat(FILE *fileIn, int id[NIDS], double *data){
  /*Reads fileIn and returns coordinates for all objects. 
    data must be allocated before.*/
  
  long i,j,Ncol;
  char line[NFIELD*NCHAR], item[NFIELD*NCHAR];
  
  i = 0;
  while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){
    if(getStrings(line,item," ",&Ncol)){
      for(j=0;j<NIDS;j++){
	data[NIDS*i+j] = getDoubleValue(item,id[j]);
      }
      i++;
    }
  }
  rewind(fileIn);
  
  return;
}

Point *dataToPoint(Config para, double *data, int N){
  /*Transforms data into array with type "Point".*/    
  long i,j,k;
  Point *result = createPoint(N,para.nboots+1);
  
  //Read file
  for(i=0;i<N;i++){
    result[i].x     = (double *)malloc(NDIM*sizeof(double));
    result[i].sinx  = (double *)malloc(2*sizeof(double));
    result[i].cosx  = (double *)malloc(2*sizeof(double));
    result[i].x[0]     = data[NIDS*i+0];
    result[i].x[1]     = data[NIDS*i+1]; 
    result[i].x[2]     = data[NIDS*i+2];
    result[i].sigz     = data[NIDS*i+3]; 
    result[i].e1       = data[NIDS*i+4]; 
    result[i].e2       = data[NIDS*i+5]; 
    result[i].w        = data[NIDS*i+6];
    result[i].cosx[0]  = cos(result[i].x[0]*PI/180.0);
    result[i].cosx[1]  = cos(result[i].x[1]*PI/180.0);
    result[i].sinx[0]  = sin(result[i].x[0]*PI/180.0);
    result[i].sinx[1]  = sin(result[i].x[1]*PI/180.0);
    //result[i].index    = i;
    result[i].boots[0] = 1;
    for(j=1;j<para.nboots+1;j++) result[i].boots[j] = 0;
    result[i].SplitDim  = 0;
  }
  
  //Bootstrap
  for(j=1;j<para.nboots+1;j++){
    for(i=0;i<N;i++){
      k = randInt(N);
      result[k].boots[j] += 1;
    }
  }
  
  return result;
}

double distAngSpher(Point *a, Point *b){
  /*Returns the distance between a and b. Spherical coordinates.*/
  double sin2_ra = 0.5*(1.0 - a->cosx[0]*b->cosx[0] - a->sinx[0]*b->sinx[0]);
  double sin2_dec = 0.5*(1.0 - a->cosx[1]*b->cosx[1] - a->sinx[1]*b->sinx[1]);
  return 2.0*asin(sqrt(MAX(EPS,sin2_dec + a->cosx[1]*b->cosx[1]*sin2_ra)))*180/PI;
}
double distAngCart(Point *a, Point *b){
  /*Returns the distance between a and b. Cartesian coordinates.*/
  return sqrt((a->x[0]-b->x[0])*(a->x[0]-b->x[0]) + (a->x[1]-b->x[1])*(a->x[1]-b->x[1]));
}

int comparePoints(const void *a,const void *b){
  /* Compares x0 coordinates of 2 points*/
 
  int SplitDim = (*(Point *)a).SplitDim;

  if ((*(Point *)a).x[SplitDim] > (*(Point *)b).x[SplitDim])
    return 1;
  else if ((*(Point *)a).x[SplitDim] < (*(Point *)b).x[SplitDim])
    return -1;
  else 
    return 0; 
}

void cpyPoints(Point *a, Point *b){
  /*Copies b into a*/
  long i;
  
  /*
  a->x     = (double *)malloc(NDIM*sizeof(double));
  a->sinx  = (double *)malloc(2*sizeof(double));
  a->cosx  = (double *)malloc(2*sizeof(double));
  */
  a->SplitDim = b->SplitDim;
  for(i=0;i<NDIM;i++){
    a->x[i] = b->x[i];
    if(i<2){
      a->sinx[i] = b->sinx[i];
      a->cosx[i] = b->cosx[i];
    }
  }
  //a->Nboots = b->Nboots;
  //a->boots  = (unsigned char *)malloc((a->Nboots)*sizeof(unsigned char));
  for(i=0;i<a->Nboots;i++){
    a->boots[i] = b->boots[i];
    }      
  //a->index = b->index;
  a->sigz  = b->sigz;
  a->e1    = b->e1;
  a->e2    = b->e2;
  a->w     = b->w;
  
}

Point *createPoint(long N, short Nboots){
  long i;
  Point *result = (Point *)malloc(N*sizeof(Point));
  for(i=0;i<N;i++){
    result[i].x      = (double *)malloc(NDIM*sizeof(double));
    result[i].sinx   = (double *)malloc(2*sizeof(double));
    result[i].cosx   = (double *)malloc(2*sizeof(double));
    result[i].Nboots = Nboots;
    result[i].boots  = (unsigned char *)malloc((Nboots)*sizeof(unsigned char));
  }
  return result;
}

void free_Point(Point *point, long N){
  long i;
  for(i=0;i<N;i++) {
    free(point[i].boots);
    free(point[i].x);
    free(point[i].sinx);
    free(point[i].cosx);
  }
  free(point);
}

void free_Node(Node *node){

  if(node->type != LEAF){
    free_Node((Node *)node->Left);
    free_Node((Node *)node->Right);
  }
  
  Node *result = (Node *)node;
  free(result->N);
  free_Point(result->point,result->Npoint);
  free(result);
  
  return;
}

int compareDoubles(const void *a,const void *b){
  /* Compares two double precision numbers*/
  if (*(double *)a > *(double *)b) 
    return 1;
  else if (*(double *)a < *(double *)b) 
    return -1;
  else 
    return 0;
}

double determineMachineEpsilon()
{
  double u, den;
  
  u = 1.0;
  do {
    u /= 2.0;
    den = 1.0 + u;
  } while(den>1.0);
  
  return(10.0 * u);
}

long determineIntError(){
  /*Equals long_max in fact.*/
  long count=0;
  count--;
  //count = (long)(-1);
  return count;
}

FILE *fopenAndCheck(const char *fileName,char *mode){
  /*Checks if fileName exists and opens it. Exits otherwise.;*/
  FILE *fileTmp = fopen(fileName,mode);
  if (fileTmp == NULL){
    printf("%s: %s not found. Exiting...\n",MYNAME,fileName);
    exit(EXIT_FAILURE);    
  }
  return fileTmp;
}

int getStrings(char *line, char *strings, char *delimit, long *N){
  /*Extract each word/number in line separated by delimit and returns 
    the array of items in strings.*/
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

int roundToNi(double a){
  /*round a to the next integer*/
  return (int)round(a);  
}
void printCount(const long *count, const long *total, int verbose){
  if(verbose){
    if((*count)%100000 == 0){
      fflush(stdout);
      printf("\b\b\b\b\b\b\b%6.2f%%",100.0*(double)(*count)/(double)(*total));
    }
  }
}

long randInt(long N){
  /*Returns an random integer between 0 and N-1*/
  return (long)((double)rand()/((double)RAND_MAX + 1)*N);
}


/*----------------------------------------------------------------*
 *Obsolete stuff                                                  *
 *----------------------------------------------------------------*/
/*
   in Npairs: ------------------------------------------
    for(i=0;i<node1->N[0];i++){
      for(j=0;j<node2->N[0];j++){
	d = para->distAng(&node1->rootPoint[node1->index[i]],&node2->rootPoint[node2->index[j]]);
	if(para->log) d = log(d);
	  k = floor((d-para->min)/para->Delta);
	  if(0 <= k && k < para->nbins){
	    for(l=0;l<para->nboots+1;l++) N[l+(para->nboots+1)*k] += node1->rootPoint[node1->index[i]].boots[l]*node2->rootPoint[node2->index[j]].boots[l];
	  }
      }
    }
    
    In createNode: ----------------------------------------
    //index of objects inside the node
    result->index = (long *)malloc(N*sizeof(long));
    for(i=0;i<N;i++){
    result->index[i] = data[i].index;
    }
  

Point meanPoints(Point *point, long N){
  //Return mean properties
  long i,j,k;
  Point result;
  double sum;
  
  //Barycenter
  for(j=0;j<NDIM;j++){
    result.x[j] = 0.0;
    for(i=0;i<N;i++) {
      result.x[j] += point[i].x[j]/(double)N;
    }
    result.cosx[j] = cos(result.x[j]*PI/180.0);
    result.sinx[j] = sin(result.x[j]*PI/180.0);
  }
  

  result.SplitDim = point[0].SplitDim;
  result.Nboots   = point[0].Nboots;
  result.boots    = (int *)malloc((result.Nboots)*sizeof(int));
  for(j=0;j<result.Nboots;j++){
    sum = 0.0;
    for(i=0;i<N;i++) {
      sum += (double)point[i].boots[j]/(double)N;
    }
    result.boots[j] = roundToNi(sum);
  }
  
  return result;
}






    */
