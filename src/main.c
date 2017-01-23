#include "main.h"

/*
 *    main.c
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2017)
 *
 *    Contributions:
 *    - the core algorithm to compute the number of pairs from a kd-tree is
 *    based on Martin Kilbinger's Athena code: http://www.cosmostat.org/software/athena/.
 *    - the galaxy-galaxy lensing core algorithm is based on Alexie Leauthaud's
 *    code (Leauthaud et al. 2010, ApJ, 709, 97).
 *    - Elinor Medezinski helped improving and correcting a number of bugs in the
 *    cross-correlation module.
 *
 *    Program to compute two-point correlation functions based on
 *    "divide and conquer" algorithms, mainly, but not limited to:
 *    - data storage in a binary tree,
 *    - approximation at large scale,
 *    - parellelization.
 *    Supports auto and cross correlations, and galaxy-galaxy lensing.
 *
 *    TODO:
 *    - check calibration factor (simply replace e by 1+m ?)
 *    - approximation at large scale for 3D and wp based on physical distance
 *    - compute tangential shear
 *    - option to take into account the East-West orientation
 *    - Do the randoms one by one or all together at once? -> Still need
 *    to decide this question. NOT URGENT
 *
 *    Known bugs:
 *
 *
 */

int main(argc,argv)
int argc;
char **argv;
{

   int rank = 0, size = 1;
   /*    MPI initialization:
    *    "rank" is the id of the current cpu in use [0: master, size-1: last cpu]
    *    "size" is the total number of cpus in use ("-np" option when calling mpirun)
    */
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   /*    global configuration -> read config file and command line
    *    note: all cpus read set para on their side because
    *    sending structures is much more complex.
    */
   Config para;
   para.size = size;
   para.rank = rank;
   initPara(argc,argv,&para);

   /*    computation time begins */
   double t0 = MPI_Wtime();

   switch (para.corr){
      case AUTO: case AUTO_WP: case AUTO_3D:
         /*    two-point autocorrelation function */
         autoCorr(para);
         break;
      case CROSS: case CROSS_WP: case CROSS_3D:
         /*    two-point cross-correlation function */
         crossCorr(para);
         break;
      case GGLENS:
         /*    galaxy-galaxy lensing two-point cross-correlation function */
         ggCorr(para);
         break;
      case NUMBER :
         /*    number counts */
         numberCount(para);
         break;
   }

   /*    computation time ends */
   MPI_Barrier(MPI_COMM_WORLD);
   if(para.verbose) fprintf(stderr, "%s: elapsed time: %5.2f s\n", MYNAME, MPI_Wtime() - t0);

   /*    end of the main program */
   MPI_Finalize();

   return EXIT_SUCCESS;
}

/*
 *    Main routines
 */

void numberCount(Config para){
   /*    Computes the number of objects
    *    for autocorrelation, each cpu correlates all nodes from the root
    *    with all nodes from a subnode (which is only a part of the tree).
    */

   int dimStart = 0;
   long i, j, k, l, n;
   double DD_sum, DR_sum, RR_sum;
   Point data, random;
   Tree dataTree, randomTree;
   Result N;
   Mask mask;
   char  fileOutName[1000];
   FILE *fileOut;

   /* read files */
   if(para.rank == MASTER){
      comment(para,"Reading fileRan1..."); random = readCat(para, para.fileRanName1, para.ran1Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random.N);}

      comment(para,"Reading fileIn1....");  data  = readCat(para, para.fileInName1, para.data1Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", data.N);}
   }

   /*    resample, build masks */
   comment(para, "Resampling...");
   Mask limits;
   if(para.rank == MASTER){

      limits.min = (double *)malloc(NDIM*sizeof(double));
      limits.max = (double *)malloc(NDIM*sizeof(double));
      setLimits(&data, &limits);
   }

   resample(&para, &random, dimStart, &mask, &limits, FIRSTCALL);

   if(para.rank == MASTER){
      free(limits.min);
      free(limits.max);
   }

   /*    send data */
   comment(para, "sending data...");
   comData(para, &random, 0, dimStart, FIRSTCALL);
   comData(para, &data  , 0, dimStart, FIRSTCALL);

   /*    grow trees */
   comment(para, "building trees...");
   dataTree   = buildTree(&para, &data, &mask, dimStart, FIRSTCALL);     freePoint(para, data);

   /*    Output tree to ascii file. format: RA DEC [w_0...w_nsamples] rank*/
   if(para.rank == 0 && (para.printTree || para.printTreeAndExit)){
      comment(para, "writing trees...");
      if(para.fits){
         sprintf(fileOutName, "!%s.data1_tree.fits", para.fileOutName);
         printTreeFits(para, fileOutName, dataTree, ROOT, 1, FIRSTCALL);
      }else{
         sprintf(fileOutName, "%s.data1_tree.ascii", para.fileOutName);
         printTree(para, fileOutName, dataTree, ROOT, 1, FIRSTCALL);
      }
   }
   comment(para, "done.\n");
   if(para.printTreeAndExit){
      return;
   }

   /*    divide and conquer */
   long nodeSlaveData = splitTree(&para, &dataTree, ROOT, para.size, FIRSTCALL);

   /*    compute number of objects */
   comment(para, "N...       "); N = Nobjects(&para, &dataTree, nodeSlaveData, FIRSTCALL);

   freeMask(para, mask);
   freeTree(para, dataTree);

   /* each slave sends the result and master sums everything up */
   comResult(para, N, para.size, 0);

   /* print out results */
   if(para.rank == MASTER){

      /*    R */
      double *R          = (double *)malloc(para.nbins*sizeof(double));
      double sum, *meanR = (double *)malloc(para.nbins*sizeof(double));
      for(i=0;i<para.nbins;i++){
         sum = N.NN[i];
         if(para.log){
            R[i] = meanR[i] = exp(para.min+para.Delta*(double)i+para.Delta/2.0);
            if(sum > 0.0) meanR[i] = exp(N.meanR[i]/sum);
         }else{
            R[i] = meanR[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
            if(sum > 0.0) meanR[i] = N.meanR[i]/sum;
         }
      }

      /*    Errors  */

      /*    mean N and errors */
      double *Nmean = (double *)malloc(para.nbins*sizeof(double));
      double *err_r = (double *)malloc(para.nbins*sizeof(double));
      double *err_p = (double *)malloc(para.nbins*sizeof(double));

      double norm, norm_N;
      switch(para.err){
         case JACKKNIFE: norm = (double)(para.nsamples - 1)/(double)(para.nsamples); break;
         case BOOTSTRAP: norm = 1.0/(double)(para.nsamples - 1); break;
         case SUBSAMPLE: norm = (double)(para.nsamples - 1); break;
      }

      for(i=0;i<para.nbins;i++){

         Nmean[i] = err_r[i] = 0.0;

         /*    1. resampling error */
         if(para.nsamples > 1){ /* mean */
            for(l=0;l<para.nsamples;l++){
               //norm_N = N.N1[0]/N.N1[l+1];
               norm_N = 1.0;
               Nmean[i] += norm_N  * N.NN[para.nbins*(l+1)+i]/(double)para.nsamples;
            }
            for(l=0;l<para.nsamples;l++){ /* dispersion */
               err_r[i] += SQUARE(Nmean[i]-norm_N*N.NN[para.nbins*(l+1)+i]);
            }
            err_r[i] = sqrt(norm*err_r[i]);
         }

         /*    2. poisson error ~1/N */
         err_p[i] = sqrt(N.NN[para.nbins*0+i]);
      }

      /*    write file out */
      fileOut = fopen(para.fileOutName, "w");
      fprintf(fileOut, "# Number counts.\n");
      switch(para.err){
         case JACKKNIFE: fprintf(fileOut, "# Resampling: jackknife (%d samples [=subvolumes]).\n", para.nsamples); break;
         case BOOTSTRAP: fprintf(fileOut, "# Resampling: bootstrap (%d subvolumes, %d samples).\n", para.nsub, para.nsamples); break;
         case SUBSAMPLE: fprintf(fileOut, "# Resampling: subsample (%d samples [=subvolumes]).\n", para.nsamples); break;
      }
      fprintf(fileOut, "#  x            N(x)    err(resamp) err(resamp-poisson) err(poisson)      <R>\n");
      for(i=0;i<para.nbins;i++){
         fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n",  R[i], N.NN[i], err_r[i], sqrt(MAX(0.0, SQUARE(err_r[i])-SQUARE(err_p[i]))), err_p[i], meanR[i]);
      }
      fclose(fileOut);

      /*    write samples */
      if(para.printSamples){
         sprintf(fileOutName, "%s.samples", para.fileOutName);
         /* TODO */
      }

      /*    write file out covariance matrix */
      if(para.cov_mat){
         double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
         sprintf(fileOutName, "%s.cov", para.fileOutName);
         fileOut = fopen(fileOutName, "w");
         for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
         for(i=0;i<para.nbins;i++){
            for(j=0;j<para.nbins;j++){
               for(l=0;l<para.nsamples;l++){
                  cov[para.nbins*i+j] += norm*(Nmean[i]-N.NN[para.nbins*(l+1)+i])
                  *(Nmean[j]-N.NN[para.nbins*(l+1)+j]);
               }
               fprintf(fileOut,"%g ", cov[para.nbins*i+j]);
            }
            fprintf(fileOut,"\n");
         }
         fclose(fileOut);
         free(cov);
      }

      free(R);
      free(meanR);
      free(Nmean);
      free(err_r);
      free(err_p);
   }

   freeResult(para, N);

   return;
}

void autoCorr(Config para){
   /*    Computes the auto correlation function.
    *    for autocorrelation, each cpu correlates all nodes from the root
    *    with all nodes from a subnode (which is only a part of the tree).
    */

   int dimStart = 0;
   long i, j, k, l, n;
   double DD_sum, DR_sum, RR_sum;
   Point data, random;
   Tree dataTree, randomTree;
   Result DD, DR, RR;
   Mask mask;
   char  fileOutName[1000];
   FILE *fileOut, *fileRR;

   /*    read files */
   if(para.rank == MASTER){

      comment(para,"Reading fileRan1..."); random = readCat(para, para.fileRanName1, para.ran1Id,  para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random.N);}

      comment(para,"Reading fileIn1....");  data  = readCat(para, para.fileInName1, para.data1Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n",data.N);}
   }

   /*    resample, build masks */
   comment(para, "Resampling...");

   Mask limits;
   if(para.rank == MASTER){
      limits.min = (double *)malloc(NDIM*sizeof(double));
      limits.max = (double *)malloc(NDIM*sizeof(double));
      setLimits(&data, &limits);
   }

   resample(&para, &random, dimStart, &mask, &limits, FIRSTCALL);

   if(para.rank == MASTER){
      free(limits.min);
      free(limits.max);
   }


   /*    send data */
   comment(para, "sending data...");
   comData(para, &random, 0, dimStart, FIRSTCALL);
   comData(para, &data  , 0, dimStart, FIRSTCALL);

   /*    grow trees */
   comment(para, "building trees...");
   randomTree = buildTree(&para, &random, &mask, dimStart, FIRSTCALL);   freePoint(para, random);
   dataTree   = buildTree(&para, &data, &mask, dimStart, FIRSTCALL);     freePoint(para, data);

   /*    divide and conquer */
   long nodeSlaveRan  = splitTree(&para, &randomTree, ROOT, para.size, FIRSTCALL);
   long nodeSlaveData = splitTree(&para, &dataTree,   ROOT, para.size, FIRSTCALL);

   /*    output tree to ascii file. format: RA DEC [w_0...w_nsamples] rank*/
   if(para.rank == 0 && (para.printTree || para.printTreeAndExit)){
      comment(para, "writing trees...");
      if(para.fits){
         sprintf(fileOutName, "!%s.data1_tree.fits", para.fileOutName);
         printTreeFits(para, fileOutName, dataTree, ROOT, 1, FIRSTCALL);
      }else{
         sprintf(fileOutName, "%s.data1_tree.ascii", para.fileOutName);
         printTree(para, fileOutName, dataTree, ROOT, 1, FIRSTCALL);
      }
   }
   comment(para, "done.\n");
   if (para.printTreeAndExit){
      return;
   }

   /*    compute pairs */
   switch (para.corr){
      case AUTO: case AUTO_3D:
      if (!strcmp(para.RRInFileName, "")){
         comment(para, "RR...       "); RR = Npairs(&para, &randomTree, ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
      }
      comment(para, "DR...       "); DR = Npairs(&para, &dataTree,   ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
      comment(para, "DD...       "); DD = Npairs(&para, &dataTree,   ROOT, &dataTree,   nodeSlaveData, FIRSTCALL);
      break;
      case AUTO_WP:
      if (!strcmp(para.RRInFileName, "")){
         comment(para, "RR(rp,pi)...       "); RR = NpairsWp(&para, &randomTree, ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
      }
      comment(para, "DR(rp,pi)...       "); DR = NpairsWp(&para, &dataTree,   ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
      comment(para, "DD(rp,pi)...       "); DD = NpairsWp(&para, &dataTree,   ROOT, &dataTree,   nodeSlaveData, FIRSTCALL);
      break;
   }

   freeMask(para, mask);
   freeTree(para, randomTree);
   freeTree(para, dataTree);

   /*    each slave sends the result and master sums everything up */
   if (!strcmp(para.RRInFileName, "")){
      comResult(para, RR, para.size, 0);
   }
   comResult(para, DR, para.size, 0);
   comResult(para, DD, para.size, 0);

   /*    RR pairs */
   if(para.rank == MASTER){
      /*       Save RR pairs */
      if (strcmp(para.RROutFileName, "")){
         fileRR = fopen(para.RROutFileName, "w");
         if(para.corr == AUTO || para.corr == AUTO_3D){
            fwrite(RR.NN, sizeof(double),    para.nbins*(para.nsamples+1), fileRR);
            fwrite(RR.N1, sizeof(double),    para.nsamples+1, fileRR);
            fwrite(RR.N2, sizeof(double),    para.nsamples+1, fileRR);
            fwrite(RR.meanR, sizeof(double), para.nbins, fileRR);
         }else if(para.corr == AUTO_WP){
            fwrite(RR.NN,    sizeof(double), para.nbins*para.nbins_pi*(para.nsamples+1), fileRR);
            fwrite(RR.N1,    sizeof(double), para.nsamples+1, fileRR);
            fwrite(RR.N2,    sizeof(double), para.nsamples+1, fileRR);
            fwrite(RR.meanR, sizeof(double), para.nbins, fileRR);
            fwrite(RR.NN_s,  sizeof(double), para.nbins*(para.nsamples+1), fileRR);
         }
         fclose(fileRR);
      }

      /*    Get RR pairs from previous run */
      /* TODO: save bins and check them at load time */
      if (strcmp(para.RRInFileName, "")){
         fileRR = fopen(para.RRInFileName, "r");
         if(para.corr == AUTO || para.corr == AUTO_3D){
            RR.NN    = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));
            RR.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
            RR.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
            RR.meanR = (double *)malloc(para.nbins*sizeof(double));

            fread(RR.NN, sizeof(double),    para.nbins*(para.nsamples+1), fileRR);
            fread(RR.N1, sizeof(double),    para.nsamples+1, fileRR);
            fread(RR.N2, sizeof(double),    para.nsamples+1, fileRR);
            fread(RR.meanR, sizeof(double), para.nbins, fileRR);
         }else if(para.corr == AUTO_WP){
            RR.NN    = (double *)malloc(para.nbins*para.nbins_pi*(para.nsamples+1)*sizeof(double));
            RR.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
            RR.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
            RR.meanR = (double *)malloc(para.nbins*sizeof(double));
            RR.NN_s  = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));

            fread(RR.NN,    sizeof(double), para.nbins*para.nbins_pi*(para.nsamples+1), fileRR);
            fread(RR.N1,    sizeof(double), para.nsamples+1, fileRR);
            fread(RR.N2,    sizeof(double), para.nsamples+1, fileRR);
            fread(RR.meanR, sizeof(double), para.nbins, fileRR);
            fread(RR.NN_s,  sizeof(double), para.nbins*(para.nsamples+1), fileRR);
         }
         fclose(fileRR);
      }
   }

   /*    print out results */
   if(para.rank == MASTER){

      /*    R */
      double *R          = (double *)malloc(para.nbins*sizeof(double));
      double sum, *meanR = (double *)malloc(para.nbins*sizeof(double));
      for(i=0;i<para.nbins;i++){
         if(para.corr == AUTO || para.corr == AUTO_3D ) sum = DD.NN[i];
         if(para.corr == AUTO_WP) {
            sum = 0.0;
            // here i -> rp, j -> pi
            for(j=0;j<para.nbins_pi;j++) sum += DD.NN[i + para.nbins*j];
         }
         if(para.log){
            R[i] = meanR[i] = exp(para.min+para.Delta*(double)i+para.Delta/2.0);
            if(sum > 0.0) meanR[i] = exp(DD.meanR[i]/sum);
         }else{
            R[i] = meanR[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
            if(sum > 0.0) meanR[i] = DD.meanR[i]/sum;
         }
      }



      /*    Errors  */
      /*    mean w(theta) and errors */
      double *wmean = (double *)malloc(para.nbins*sizeof(double));
      double *err_r = (double *)malloc(para.nbins*sizeof(double));
      double *err_p = (double *)malloc(para.nbins*sizeof(double));

      double norm;
      switch(para.err){
         case JACKKNIFE: norm = (double)(para.nsamples - 1)/(double)(para.nsamples); break;
         case BOOTSTRAP: norm = 1.0/(double)(para.nsamples - 1); break;
      }

      for(i=0;i<para.nbins;i++){
         wmean[i] = err_r[i] = 0.0;
         /*    1. resampling error */
         if(para.nsamples > 1){ /*  mean */
            for(l=0;l<para.nsamples;l++){
               switch (para.corr){
                  case AUTO: case AUTO_3D:
                     wmean[i] += wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1)/(double)para.nsamples;
                     break;
                  case AUTO_WP:
                     wmean[i] += wp(para, para.estimator, DD, RR, DR, DR, -1, i, l+1)/(double)para.nsamples;
                     // wmean[i] += xis(para, para.estimator, DD, RR, DR, DR, i, l+1)/(double)para.nsamples;
                     break;
               }
            }
            for(l=0;l<para.nsamples;l++){ /* dispersion */
               switch (para.corr){
                  case AUTO: case AUTO_3D:
                     err_r[i] += SQUARE(wmean[i]-wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1));
                     break;
                  case AUTO_WP:
                     err_r[i] += SQUARE(wmean[i]-wp(para, para.estimator, DD, RR, DR, DR, -1, i, l+1));
                     //err_r[i] += SQUARE(wmean[i]-xis(para, para.estimator, DD, RR, DR, DR, i, l+1));
                  break;
               }
            }
            err_r[i] = sqrt(norm*err_r[i]);
         }
         /*    2. poisson error ~1/N */
         switch (para.corr){
            case AUTO:  case AUTO_3D:
               if(DD.NN[para.nbins*0+i] > 0 && RR.NN[para.nbins*0+i] > 0){
                  err_p[i] = ABS(1.0+wTheta(para, para.estimator, DD, RR, DR, DR, i, 0))*(1.0/sqrt((double)DD.NN[para.nbins*0+i]) + 1.0/sqrt((double)RR.NN[para.nbins*0+i]));
               }else{
                  err_p[i] = 0.0;
               }
               break;
            case AUTO_WP:
               // i,j are the bin indexes. here i : rp, j : pi
               DD_sum = 0;
               RR_sum = 0;
               j      = 0;
               // while(R[j] < para.pi_max && j < para.nbins_pi){ /* sum over pi integration */
               for(j=0;j<para.nbins_pi;j++){
                  DD_sum += DD.NN[i + para.nbins*j];
                  RR_sum += RR.NN[i + para.nbins*j];
                  //j++;
               }
               if(DD_sum > 0 && RR_sum > 0){
                  err_p[i] = ABS(1.0+wp(para, para.estimator, DD, RR, DR, DR, -1, i, 0))*(1.0/sqrt((double)DD_sum + 1.0/sqrt((double)RR_sum)));
               }else{
                  err_p[i] = 0.0;
               }
               break;
         }
      }

      /*    write file out */
      fileOut = fopen(para.fileOutName, "w");
      switch(para.estimator){
         case LS:      fprintf(fileOut, "# Auto-correlation, Landy & Szalay estimator.\n"); break;
         case NAT:     fprintf(fileOut, "# Auto-correlation, Natural estimator.\n");        break;
         case HAM:     fprintf(fileOut, "# Auto-correlation, Hamilton estimator.\n");       break;
         case PEEBLES: fprintf(fileOut, "# Auto-correlation, Peebles estimator.\n");        break;
      }
      switch(para.err){
         case JACKKNIFE: fprintf(fileOut, "# Resampling: jackknife (%d samples [=subvolumes]).\n", para.nsamples); break;
         case BOOTSTRAP: fprintf(fileOut, "# Resampling: bootstrap (%d subvolumes, %d samples).\n", para.nsub, para.nsamples); break;
      }
      switch(para.corr){
         /* If auto correlation  */
         case AUTO:  case AUTO_3D:
            fprintf(fileOut, "#  R            w       err(resamp) err(resamp-poisson) err(poisson)      <R>                DD                DR                RR      Ndata             Nrandom\n");
            /*
            if(para.proj == COMO){
               fprintf(fileOut, "#  R(comoving)  w       err(resamp) err(resamp-poisson) err(poisson)      <R>                DD                DR                RR      Ndata             Nrandom\n");
            }else if(para.proj == PHYS){
               fprintf(fileOut, "#  R(physical)  w       err(resamp) err(resamp-poisson) err(poisson)      <R>                DD                DR                RR      Ndata             Nrandom\n");
            }else{
               fprintf(fileOut, "#  theta        w       err(resamp) err(resamp-poisson) err(poisson)      <R>                DD                DR                RR      Ndata             Nrandom\n");
            }
            */
            for(i=0;i<para.nbins;i++){
               fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %17.5f %17.5f %17.5f %17.5f %17.5f\n",
                  R[i], wTheta(para, para.estimator, DD, RR, DR, DR, i, 0),
                  err_r[i], sqrt(MAX(0.0, SQUARE(err_r[i])-SQUARE(err_p[i]))), err_p[i], meanR[i],
                  DD.NN[i], DR.NN[i], RR.NN[i],  DD.N1[0],  RR.N1[0]);
            }
            break;

         /* If wp(rp) auto correlation  */
         case AUTO_WP:
            fprintf(fileOut, "# pi upper limit integration: %5.2f Mpc.\n", para.pi_max);
            fprintf(fileOut, "# Attention: sum(...) are the integrated sums of pairs along pi,\n");
            fprintf(fileOut, "# given for reference ONLY. No combination of these would lead to wp.\n");
            if(para.proj == COMO){
               fprintf(fileOut, "# rp (comoving) wp    err(resamp) err(resamp-poisson) err(poisson)        <R>           sum(DD)           sum(DR)           sum(RR)             Ndata           Nrandom\n");
            }else if(para.proj == PHYS){
               fprintf(fileOut, "# rp (physical) wp    err(resamp) err(resamp-poisson) err(poisson)        <R>           sum(DD)           sum(DR)           sum(RR)             Ndata           Nrandom\n");
            }
            for(i=0;i<para.nbins;i++){
               RR_sum = 0;
               DR_sum = 0;
               DD_sum = 0;
               j      = 0;
               // while(R[j] < para.pi_max && j < para.nbins_pi){  /* sum over pi integration - ONLY informative */
               for(j=0;j<para.nbins_pi;j++){
                  RR_sum += RR.NN[i + para.nbins*j];
                  DR_sum += DR.NN[i + para.nbins*j];
                  DD_sum += DD.NN[i + para.nbins*j];
                  //j++;
               }

               /* TEST: xis(para, para.estimator, DD, RR, DR, DR, i, 0),*/

               fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %17.5f %17.5f %17.5f %17.5f %17.5f\n",
                  R[i], wp(para, para.estimator, DD, RR, DR, DR, -1, i, 0),
                  err_r[i], sqrt(MAX(0.0, SQUARE(err_r[i])-SQUARE(err_p[i]))),  err_p[i], meanR[i],
                  DD_sum, DR_sum, RR_sum,  DD.N1[0],  RR.N1[0]);
            }
            break;
      }
      fclose(fileOut);

      /*    write file out xi(rp,pi) */
      if(para.corr == AUTO_WP && para.xi){
         sprintf(fileOutName, "%s.xi", para.fileOutName);
         fileOut = fopen(fileOutName,"w");
         fprintf(fileOut, "# Rows: rp, columns: pi\n#");
         for(i=0;i<para.nbins;i++) fprintf(fileOut, "%12.7f ", R[i]);
         fprintf(fileOut, "\n");
         for(i=0;i<para.nbins;i++){
            for(j=0;j<para.nbins_pi;j++){
               fprintf(fileOut, "%12.7f ", wp(para, para.estimator, DD, RR, DR, DR, j, i, 0));
            }
            fprintf(fileOut, "\n");
         }
         fclose(fileOut);
      }

      /*    write samples */
      if(para.printSamples){
         sprintf(fileOutName, "%s.samples", para.fileOutName);
         fileOut = fopen(fileOutName, "w");

         /*    Print resampling method */
         switch(para.err){
            case JACKKNIFE: fprintf(fileOut, "# Resampling: jackknife (%d samples [=subvolumes]).\n", para.nsamples); break;
            case BOOTSTRAP: fprintf(fileOut, "# Resampling: bootstrap (%d subvolumes, %d samples).\n", para.nsub, para.nsamples); break;
         }

         /*    Print bin values */
         fprintf(fileOut, "#");
         for(i=0;i<para.nbins;i++){
            fprintf( fileOut, "%12.7f ", R[i]);
         }
         fprintf(fileOut, "\n");

         switch(para.corr){
            case AUTO: case AUTO_3D:
               for(l=0;l<para.nsamples;l++){
                  for(i=0;i<para.nbins;i++){
                     fprintf(fileOut, "%12.7f ", wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1) );
                  }
                  fprintf(fileOut, "\n");
               }
               break;
            case AUTO_WP:
               for(l=0;l<para.nsamples;l++){
                  for(i=0;i<para.nbins;i++){
                     fprintf(fileOut, "%12.7f ", wp(para, para.estimator, DD, RR, DR, DR, -1, i, l+1) );
                  }
                  fprintf(fileOut, "\n");
               }
               break;
         }
         fclose(fileOut);
      }


      /*    write file out covariance matrix */
      if(para.cov_mat){
         double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
         for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
         sprintf(fileOutName, "%s.cov", para.fileOutName);
         fileOut = fopen(fileOutName, "w");
         for(i=0;i<para.nbins;i++){
            for(j=0;j<para.nbins;j++){
               switch(para.corr){
                  case AUTO:  case AUTO_3D:
                     for(l=0;l<para.nsamples;l++){
                        cov[para.nbins*i+j] += norm*(wmean[i]-wTheta(para, para.estimator, DD, RR, DR, DR, i, l+1))
                        *(wmean[j]-wTheta(para, para.estimator, DD, RR, DR, DR, j, l+1));
                     }
                     break;
                  case AUTO_WP:
                     for(l=0;l<para.nsamples;l++){
                        cov[para.nbins*i+j] += norm*(wmean[i]-wp(para, para.estimator, DD, RR, DR, DR, -1, i, l+1))
                        *(wmean[j]-wp(para, para.estimator, DD, RR, DR, DR, -1, j, l+1));
                     }
                     break;
               }
               fprintf(fileOut,"%g ", cov[para.nbins*i+j]);
            }
            fprintf(fileOut,"\n");
         }
         fclose(fileOut);
         free(cov);
      }

      free(R);
      free(meanR);
      free(wmean);
      free(err_r);
      free(err_p);
   }

   if(!strcmp(para.RRInFileName, "") || para.rank == MASTER){
      freeResult(para, RR);
   }
   freeResult(para, DR);
   freeResult(para, DD);

   return;
}

void crossCorr(Config para){
   /*    Computes the cross-correlation function.
   *     for autocorrelation, each cpu correlates all nodes from the root
   *     with all nodes from a subnode (which is only a part of the tree).
   */

   int swapped = 0, dimStart = 0;
   long i, j, k, l, n;
   double  D1D2_sum, D1R2_sum, D2R1_sum, R1R2_sum;
   Point data1, data2, random1, random2, tmp;
   Tree dataTree1, dataTree2, randomTree1, randomTree2;
   Result D1D2, D1R2, D2R1, R1R2;
   Mask mask;
   char  fileOutName[1000];
   FILE *fileOut, *fileRR;

   /*    read files */
   if(para.rank == MASTER){
      comment(para,"Reading fileRan1..."); random1 = readCat(para, para.fileRanName1, para.ran1Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random1.N);}

      comment(para,"Reading fileRan2..."); random2 = readCat(para, para.fileRanName2, para.ran2Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n", random2.N);}

      comment(para,"Reading fileIn1..."); data1 = readCat(para, para.fileInName1, para.data1Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n",data1.N);}

      comment(para,"Reading fileIn2..."); data2 = readCat(para, para.fileInName2, para.data2Id, para.weighted);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd objects found).\n",data2.N);}

      /* swap "1" and "2" if "1" is larger (to improve memory management)  */
      /*
      if(random1.N > random2.N){
         comment(para, "ATTENTION: fileRan1 is larger than fileRan2. Swapping \"1\" and \"2\" to save memory...\n");
         tmp = data1;   data1   = data2;   data2   = tmp;
         tmp = random1; random1 = random2; random2 = tmp;
         swapped = 1;
      }
      */
   }

   /*    resample, build mask from random file 1 */
   comment(para, "Resampling...");

   Mask limits;

   if(para.rank == MASTER){

      limits.min = (double *)malloc(NDIM*sizeof(double));
      limits.max = (double *)malloc(NDIM*sizeof(double));
      setLimits(&data1, &limits);
   }

   resample(&para, &random1, dimStart, &mask, &limits, FIRSTCALL);

   if(para.rank == MASTER){
      free(limits.min);
      free(limits.max);
   }

   /*    send data */
   comment(para, "sending data...");
   comData(para, &random1, 0, dimStart, FIRSTCALL);
   //comData(para, &random2, para.size, dimStart, FIRSTCALL); /* 2 is partitioned among cpus */
   comData(para, &random2, 0, dimStart, FIRSTCALL);
   comData(para, &data1, 0, dimStart, FIRSTCALL);
   comData(para, &data2, 0, dimStart, FIRSTCALL);


   /*    grow trees */
   comment(para, "building trees...");
   randomTree1 = buildTree(&para, &random1, &mask, dimStart, FIRSTCALL);   freePoint(para, random1);
   randomTree2 = buildTree(&para, &random2, &mask, dimStart, FIRSTCALL);   freePoint(para, random2);
   dataTree1   = buildTree(&para, &data1, &mask, dimStart, FIRSTCALL);     freePoint(para, data1);
   dataTree2   = buildTree(&para, &data2, &mask, dimStart, FIRSTCALL);     freePoint(para, data2);

   /*    Ouput tree to ascii file. format: RA DEC [w_0...w_nsamples] rank */
   if(para.rank == 0 && (para.printTree || para.printTreeAndExit)){
      comment(para, "writing trees...");
      if(para.fits){
         sprintf(fileOutName, "!%s.data1_tree.fits", para.fileOutName);
         printTreeFits(para, fileOutName, dataTree1, ROOT, 1, FIRSTCALL);
      }else{
         sprintf(fileOutName, "%s.data1_tree.ascii", para.fileOutName);
         printTree(para, fileOutName, dataTree1, ROOT, 1, FIRSTCALL);
      }
      if(para.fits){
         sprintf(fileOutName, "!%s.data2_tree.fits", para.fileOutName);
         printTreeFits(para, fileOutName, dataTree2, ROOT, 1, FIRSTCALL);
      }else{
         sprintf(fileOutName, "%s.data2_tree.ascii", para.fileOutName);
         printTree(para, fileOutName, dataTree2, ROOT, 1, FIRSTCALL);
      }
   }
   comment(para, "done.\n");
   if(para.printTreeAndExit){
      return;
   }

   /*    divide and conquer */

   long nodeSlaveRan1  = splitTree(&para, &randomTree1, ROOT, para.size, FIRSTCALL);
   long nodeSlaveRan2  = splitTree(&para, &randomTree2, ROOT, para.size, FIRSTCALL);
   // long nodeSlaveData1 = splitTree(&para, &dataTree1,   ROOT, para.size, FIRSTCALL);
   long nodeSlaveData2 = splitTree(&para, &dataTree2,   ROOT, para.size, FIRSTCALL);

   /* end of the main program */
   //MPI_Finalize();
   //exit(-1);

   /*    compute pairs */
   switch (para.corr){
      case CROSS: case CROSS_3D:
         if (!strcmp(para.RRInFileName, "")){
            comment(para, "R1R2...       "); R1R2 = Npairs(&para, &randomTree1, ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         }

         // DEBUGGING - which one is correct for Peebles? D1R2 or D1R2?
         // comment(para, "D1R2...       "); D1R2 = Npairs(&para, &dataTree1,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
         // comment(para, "D2R1...       "); D2R1 = Npairs(&para, &dataTree2,   ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         comment(para, "D1R2...       "); D1R2 = Npairs(&para, &dataTree1,   ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         comment(para, "D2R1...       "); D2R1 = Npairs(&para, &dataTree2,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
         comment(para, "D1D2...       "); D1D2 = Npairs(&para, &dataTree1,   ROOT, &dataTree2,   nodeSlaveData2, FIRSTCALL);
         break;
      case CROSS_WP:
         if (!strcmp(para.RRInFileName, "")){
            comment(para, "R1R2(rp,pi)...       "); R1R2 = NpairsWp(&para, &randomTree1, ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         }

         // DEBUGGING - which one is correct for Peebles? D1R2 or D1R2?
         // comment(para, "D1R2(rp,pi)...       "); D1R2 = NpairsWp(&para, &dataTree1,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
         // comment(para, "D2R1(rp,pi)...       "); D2R1 = NpairsWp(&para, &dataTree2,   ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         comment(para, "D1R2(rp,pi)...       "); D1R2 = NpairsWp(&para, &dataTree1,   ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         comment(para, "D2R1(rp,pi)...       "); D2R1 = NpairsWp(&para, &dataTree2,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
         comment(para, "D1D2(rp,pi)...       "); D1D2 = NpairsWp(&para, &dataTree1,   ROOT, &dataTree2,   nodeSlaveData2, FIRSTCALL);
         break;
   }

   freeMask(para, mask);
   freeTree(para, randomTree1);
   freeTree(para, randomTree2);
   freeTree(para, dataTree1);
   freeTree(para, dataTree2);

   /*    each slave sends the result and master sums everything up */
   if (!strcmp(para.RRInFileName, "")){
      comResult(para, R1R2, para.size, 0); /*   "1" to tell MASTER to sum up the total number         */
   }
   comResult(para, D1R2, para.size, 0); /*   of objects since R2 has been partitionned among cpus  */
   comResult(para, D2R1, para.size, 0);
   comResult(para, D1D2, para.size, 0);

   /*    RR pairs */
   if(para.rank == MASTER){
      /*    Save RR pairs */
      if (strcmp(para.RROutFileName, "")){
         fileRR = fopen(para.RROutFileName, "w");
         if(para.corr == CROSS || para.corr == CROSS_3D){
            fwrite(R1R2.NN, sizeof(double),    para.nbins*(para.nsamples+1), fileRR);
            fwrite(R1R2.N1, sizeof(double),    para.nsamples+1, fileRR);
            fwrite(R1R2.N2, sizeof(double),    para.nsamples+1, fileRR);
            fwrite(R1R2.meanR, sizeof(double), para.nbins, fileRR);
         }else if(para.corr == CROSS_WP){
            fwrite(R1R2.NN,    sizeof(double), para.nbins*para.nbins_pi*(para.nsamples+1), fileRR);
            fwrite(R1R2.N1,    sizeof(double), para.nsamples+1, fileRR);
            fwrite(R1R2.N2,    sizeof(double), para.nsamples+1, fileRR);
            fwrite(R1R2.meanR, sizeof(double), para.nbins, fileRR);
            fwrite(R1R2.NN_s,  sizeof(double), para.nbins*(para.nsamples+1), fileRR);
         }
         fclose(fileRR);
      }

      /*    Get RR pairs from previous run */
      if (strcmp(para.RRInFileName, "")){
         fileRR = fopen(para.RRInFileName, "r");
         if(para.corr == CROSS || para.corr == CROSS_3D){
            R1R2.NN    = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));
            R1R2.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
            R1R2.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
            R1R2.meanR = (double *)malloc(para.nbins*sizeof(double));

            fread(R1R2.NN, sizeof(double),    para.nbins*(para.nsamples+1), fileRR);
            fread(R1R2.N1, sizeof(double),    para.nsamples+1, fileRR);
            fread(R1R2.N2, sizeof(double),    para.nsamples+1, fileRR);
            fread(R1R2.meanR, sizeof(double), para.nbins, fileRR);
         }else if(para.corr == CROSS_WP){
            R1R2.NN    = (double *)malloc(para.nbins*para.nbins_pi*(para.nsamples+1)*sizeof(double));
            R1R2.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
            R1R2.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
            R1R2.meanR = (double *)malloc(para.nbins*sizeof(double));
            R1R2.NN_s  = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));

            fread(R1R2.NN,    sizeof(double), para.nbins*para.nbins_pi*(para.nsamples+1), fileRR);
            fread(R1R2.N1,    sizeof(double), para.nsamples+1, fileRR);
            fread(R1R2.N2,    sizeof(double), para.nsamples+1, fileRR);
            fread(R1R2.meanR, sizeof(double), para.nbins, fileRR);
            fread(R1R2.NN_s,  sizeof(double), para.nbins*(para.nsamples+1), fileRR);
         }
         fclose(fileRR);
      }
   }

   /*    print out results */
   if(para.rank == MASTER){

      /* R */
      double *R          = (double *)malloc(para.nbins*sizeof(double));
      double sum, *meanR = (double *)malloc(para.nbins*sizeof(double));
      for(i=0;i<para.nbins;i++){
         if(para.corr == CROSS || para.corr == CROSS_3D) sum = (double)D1D2.NN[i];
         if(para.corr == CROSS_WP) {
            sum = 0.0;
            for(j=0;j<para.nbins_pi;j++) sum += D1D2.NN[i + para.nbins*j];
         }
         if(para.log){
            R[i] = meanR[i] = exp(para.min+para.Delta*(double)i+para.Delta/2.0);
            if(sum > 0.0) meanR[i] = exp(D1D2.meanR[i]/sum);
         }else{
            R[i] = meanR[i] = para.min+para.Delta*(double)i+para.Delta/2.0;
            if(sum > 0.0) meanR[i] = D1D2.meanR[i]/sum;
         }
      }

      /*    mean w(theta) and errors */
      double *wmean = (double *)malloc(para.nbins*sizeof(double));
      double *err_r = (double *)malloc(para.nbins*sizeof(double));
      double *err_p = (double *)malloc(para.nbins*sizeof(double));

      double norm;
      switch(para.err){
         case JACKKNIFE: norm = (double)(para.nsamples - 1)/(double)(para.nsamples); break;
         case BOOTSTRAP: norm = 1.0/(double)(para.nsamples - 1); break;
      }

      /*    Errors */
      for(i=0;i<para.nbins;i++){
         wmean[i] = err_r[i] = 0.0;
         /*    1. resampling error */
         if(para.nsamples > 1){
            for(l=0;l<para.nsamples;l++){ /*  mean */
               switch (para.corr){
                  case CROSS: case CROSS_3D:
                     wmean[i] += wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, i, l+1)/(double)para.nsamples;
                     break;
                  case CROSS_WP:
                     wmean[i] += wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, i, l+1)/(double)para.nsamples;
                     break;
               }
            }
            for(l=0;l<para.nsamples;l++){ /*  dispersion */
               switch (para.corr){
                  case CROSS: case CROSS_3D:
                     err_r[i] += SQUARE(wmean[i]-wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, i, l+1));
                     break;
                  case CROSS_WP:
                     err_r[i] += SQUARE(wmean[i]-wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, i, l+1));
                     break;
               }
            }
            err_r[i] = sqrt(norm*err_r[i]);
         }
         /* 2.    poisson error ~1/N */
         switch (para.corr){
            case CROSS: case CROSS_3D:
               if(D1D2.NN[para.nbins*0+i] > 0 && R1R2.NN[para.nbins*0+i] > 0){
                  err_p[i] = ABS(1.0+wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, i, 0))*(1.0/sqrt((double)D1D2.NN[para.nbins*0+i]) + 1.0/sqrt((double)R1R2.NN[para.nbins*0+i]));
               }else{
                  err_p[i] = 0.0;
               }
               break;
            case CROSS_WP:
               D1D2_sum = 0;
               R1R2_sum = 0;
               j        = 0;
               // while(R[j] < para.pi_max && j < para.nbins){ /* sum over pi integration */
               for(j=0;j<para.nbins_pi;j++){
                  D1D2_sum += D1D2.NN[i + para.nbins*j];
                  R1R2_sum += R1R2.NN[i + para.nbins*j];
                  // j++;
               }
               if(D1D2_sum > 0 && R1R2_sum > 0){
                  err_p[i] = ABS(1.0+wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, i, 0))*(1.0/sqrt((double)D1D2_sum + 1.0/sqrt((double)R1R2_sum)));
               }else{
                  err_p[i] = 0.0;
               }
               break;
         }
      }

      /*    write file out */
      fileOut = fopen(para.fileOutName,"w");
      if(swapped) fprintf(fileOut, "# ATTENTION: \"1\" and \"2\" have been swapped to save memory.\n");
      switch(para.estimator){
         case LS:      fprintf(fileOut, "# Cross-correlation. Landy & Szalay estimator\n"); break;
         case NAT:     fprintf(fileOut, "# Cross-correlation. Natural estimator\n");        break;
         case HAM:     fprintf(fileOut, "# Cross-correlation. Hamilton estimator\n");       break;
         case PEEBLES: fprintf(fileOut, "# Cross-correlation, Peebles estimator.\n");       break;
      }
      switch(para.err){
         case JACKKNIFE: fprintf(fileOut, "# Resampling: jackknife (%d samples [=subvolumes]).\n", para.nsamples); break;
         case BOOTSTRAP: fprintf(fileOut, "# Resampling: bootstrap (%d subvolumes, %d samples).\n", para.nsub, para.nsamples); break;
      }
      switch(para.corr){
         /*    If cross correlation */
         case CROSS: case CROSS_3D:
            fprintf(fileOut, "#  R            w       err(resamp) err(resamp-poisson) err(poisson)               <R>              D1D2             D1R2             D2R1          R1R2       Ndata1     Nrandom1   Ndata2     Nrandom2\n");

            /*
            if(para.proj == COMO){
               fprintf(fileOut, "#  R(comoving)  w       err(resamp) err(resamp-poisson) err(poisson)               <R>              D1D2             D1R2             D2R1          R1R2       Ndata1     Nrandom1   Ndata2     Nrandom2\n");
            }else if(para.proj == PHYS){
               fprintf(fileOut, "#  R(physical)  w       err(resamp) err(resamp-poisson) err(poisson)               <R>              D1D2             D1R2             D2R1          R1R2       Ndata1     Nrandom1   Ndata2     Nrandom2\n");
            }else{
               fprintf(fileOut, "#  theta        w       err(resamp) err(resamp-poisson) err(poisson)               <R>              D1D2             D1R2             D2R1          R1R2       Ndata1     Nrandom1   Ndata2     Nrandom2\n");
            }
            */
            for(i=0;i<para.nbins;i++){
               fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %17.5f %17.5f %17.5f %17.5f %17.5f %17.5f %17.5f %17.5f\n",
               R[i], wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, i, 0),
               err_r[i], sqrt(MAX(0.0, SQUARE(err_r[i])-SQUARE(err_p[i]))), err_p[i], meanR[i],
               D1D2.NN[i], D1R2.NN[i], D2R1.NN[i], R1R2.NN[i],  D1D2.N1[0],  R1R2.N1[0],  D1D2.N2[0],  R1R2.N2[0]);
            }
            break;
         /*    If wp(rp) cross correlation */
         case CROSS_WP:
            fprintf(fileOut, "# pi upper limit integration: %5.2f Mpc.\n", para.pi_max);
            fprintf(fileOut, "# Attention: sum(...) are the integrated sums of pairs along pi,\n");
            fprintf(fileOut, "# given for reference ONLY. No combination of these would lead to wp.\n");
            if(para.proj == COMO){
               fprintf(fileOut, "#  rp(comoving) wp      err(resamp) err(resamp-poisson) err(poisson)     <rp>         sum(D1D2)         sum(D1R2)         sum(D2R1)         sum(R1R2)            Ndata1          Nrandom1            Ndata2          Nrandom2\n");
            }else if(para.proj == PHYS){
               fprintf(fileOut, "#  rp(physical) wp      err(resamp) err(resamp-poisson) err(poisson)     <rp>         sum(D1D2)         sum(D1R2)         sum(D2R1)         sum(R1R2)            Ndata1          Nrandom1            Ndata2          Nrandom2\n");
            }

            for(i=0;i<para.nbins;i++){
               R1R2_sum = 0;
               D1R2_sum = 0;
               D2R1_sum = 0;
               D1D2_sum = 0;
               j        = 0;
               // while(R[j] < para.pi_max && j < para.nbins){  /* sum over pi integration */
               for(j=0;j<para.nbins_pi;j++){
                  R1R2_sum += R1R2.NN[i + para.nbins*j];
                  D1R2_sum += D1R2.NN[i + para.nbins*j];
                  D2R1_sum += D2R1.NN[i + para.nbins*j];
                  D1D2_sum += D1D2.NN[i + para.nbins*j];
                  // j++;
               }
               fprintf(fileOut, "%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f %17.5f %17.5f %17.5f %17.5f %17.5f %17.5f %17.5f %17.5f\n",
               R[i], wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, i, 0),
               err_r[i], sqrt(MAX(0.0, SQUARE(err_r[i])-SQUARE(err_p[i]))), err_p[i], meanR[i],
               D1D2_sum, D1R2_sum, D2R1_sum, R1R2_sum,  D1D2.N1[0], R1R2.N1[0],  D1D2.N2[0],  R1R2.N2[0]);
            }
            break;
      }
      fclose(fileOut);

      /*    write file out xi(rp,pi) */
      if(para.corr == CROSS_WP && para.xi){
         sprintf(fileOutName, "%s.xi", para.fileOutName);
         fileOut = fopen(fileOutName,"w");
         fprintf(fileOut, "# Rows: rp, columns: pi\n#");
         for(i=0;i<para.nbins;i++) fprintf(fileOut, "%12.7f ", R[i]);
         fprintf(fileOut, "\n");
         for(i=0;i<para.nbins;i++){
            for(j=0;j<para.nbins_pi;j++){
               fprintf(fileOut, "%12.7f ", wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, j, i, 0));
            }
            fprintf(fileOut, "\n");
         }
         fclose(fileOut);
      }


      /*    write samples */
      if(para.printSamples){
         sprintf(fileOutName, "%s.samples", para.fileOutName);
         fileOut = fopen(fileOutName, "w");

         /* Print resampling method */
         switch(para.err){
            case JACKKNIFE: fprintf(fileOut, "# Resampling: jackknife (%d samples [=subvolumes]).\n", para.nsamples); break;
            case BOOTSTRAP: fprintf(fileOut, "# Resampling: bootstrap (%d subvolumes, %d samples).\n", para.nsub, para.nsamples); break;
         }

         /* Print bin values */
         fprintf(fileOut, "#");
         for(i=0;i<para.nbins;i++){
            fprintf( fileOut, "%12.7f ", R[i]);
         }
         fprintf(fileOut, "\n");

         switch(para.corr){
            case CROSS: case CROSS_3D:
               for(l=0;l<para.nsamples;l++){
                  for(i=0;i<para.nbins;i++){
                     fprintf(fileOut, "%12.7f ", wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, i, l+1) );
                  }
                  fprintf(fileOut, "\n");
               }
               break;
            case CROSS_WP:
               for(l=0;l<para.nsamples;l++){
                  for(i=0;i<para.nbins;i++){
                     fprintf(fileOut, "%12.7f ", wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, i, l+1) );
                  }
                  fprintf(fileOut, "\n");
               }
               break;
         }

         fclose(fileOut);
      }

      /*    write file out covariance matrix */
      if(para.cov_mat){
         double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
         for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
         sprintf(fileOutName, "%s.cov", para.fileOutName);
         fileOut = fopen(fileOutName, "w");
         for(i=0;i<para.nbins;i++){
            for(j=0;j<para.nbins;j++){
               switch(para.corr){
                  case CROSS: case CROSS_3D:
                     for(l=0;l<para.nsamples;l++){
                        cov[para.nbins*i+j] += norm*(wmean[i]-wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, i, l+1))
                        *(wmean[j]-wTheta(para, para.estimator, D1D2, R1R2, D1R2, D2R1, j, l+1));
                     }
                     break;
                  case CROSS_WP:
                     for(l=0;l<para.nsamples;l++){
                        cov[para.nbins*i+j] += norm*(wmean[i]-wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, i, l+1))
                        *(wmean[j]-wp(para, para.estimator, D1D2, R1R2, D1R2, D2R1, -1, j, l+1));
                     }
                     break;
               }
               fprintf(fileOut,"%g ", cov[para.nbins*i+j]);
            }
            fprintf(fileOut,"\n");
         }
         fclose(fileOut);
         free(cov);
      }

      free(R);
      free(wmean);
      free(err_r);
      free(err_p);
   }

   if(!strcmp(para.RRInFileName, "") || para.rank == MASTER){
      freeResult(para, R1R2);
   }
   freeResult(para, D1R2);
   freeResult(para, D2R1);
   freeResult(para, D1D2);

   return;
}

void ggCorr(Config para){
   /*    Computes the auto correlation function.
    *    for autocorrelation, each cpu correlates all nodes from the root
    *    with all nodes from a subnode (which is only a part of the tree).
    */

   int dimStart = 0;
   long i, j, k, l, n;
   Point lens, source;
   Tree lensTree, sourceTree;
   Mask mask;
   Result result;
   char  fileOutName[1000];
   FILE *fileOut;

   /*    read files */
   if(para.rank == MASTER){
      comment(para,"Reading sources file..."); source = readCat(para, para.fileInName2, para.data2Id, 0);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd sources found).\n", source.N);}

      comment(para,"Reading lenses file...."); lens   = readCat(para, para.fileInName1, para.data1Id, 0);
      if(para.verbose){fflush(stderr); fprintf(stderr,"(%zd lenses found).\n", lens.N);}
   }

   /*  resampling = build masks TODO: check that */
   comment(para, "Resampling...");
   Mask limits;
   if(para.rank == MASTER){
      limits.min = (double *)malloc(NDIM*sizeof(double));
      limits.max = (double *)malloc(NDIM*sizeof(double));
      setLimits(&source, &limits);
   }


   resample(&para, &source, dimStart, &mask, &limits, FIRSTCALL);

   if(para.rank == MASTER){
      free(limits.min);
      free(limits.max);
   }

   /*    send data. Source catalogue is partitioned among cpus */
   comment(para, "sending data...");

   comData(para, &source, para.size, dimStart, FIRSTCALL);
   comData(para, &lens  , 0, dimStart, FIRSTCALL);

   /*    grow trees */
   comment(para, "building trees...");
   sourceTree = buildTree(&para, &source, &mask, dimStart, FIRSTCALL);   freePoint(para, source);
   lensTree   = buildTree(&para, &lens,   &mask, dimStart, FIRSTCALL);   freePoint(para, lens);

   /*    TODO */
   if(para.rank == 0 && (para.printTree || para.printTreeAndExit)){
      if(para.verbose) fprintf(stderr,"\n%s: **WARNING** printTree not supported for gglens.\n", MYNAME);
   }
   comment(para, "done.\n");
   if(para.printTreeAndExit){
      return;
   }

   comment(para, "Correlating lenses with sources...       ");
   result = gg(&para, &lensTree, ROOT, &sourceTree, ROOT, FIRSTCALL);

   freeMask(para, mask);
   freeTree(para, sourceTree);
   freeTree(para, lensTree);

   /*    each slave sends the result and master sums everything up */
   comResult(para, result, para.size, 0);

   /*    print out results */
   if(para.rank == MASTER){

      /*    resampling errors */
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

      double sign, R, meanR;
      fileOut = fopen(para.fileOutName,"w");
      fprintf(fileOut, "# Gal-gal lensing. Delta_Sigma(R) vs R, linear approximation\n");
      switch(para.err){
         case JACKKNIFE: fprintf(fileOut, "# Resampling: jackknife (%d samples [=subvolumes]).\n", para.nsamples); break;
         case BOOTSTRAP: fprintf(fileOut, "# Resampling: bootstrap (%d subvolumes, %d samples).\n", para.nsub, para.nsamples); break;
      }
      switch(para.proj){
         case PHYS:  fprintf(fileOut, "# Coordinates: physical\n"); break;
         case COMO:  fprintf(fileOut, "# Coordinates: comoving\n"); break;
      }
      fprintf(fileOut, "# Cosmolgy: H0 = %g, Omega_M = %g, Omega_L = %g\n", para.a[0], para.a[1], para.a[2]);
      if(para.calib){
         sign = 1.0;
         fprintf(fileOut, "#  R(Mpc)  calib_factor    err(weights) err(resampling) Nsources       <R>          e2\n");
      }else{
         fprintf(fileOut, "#  R(Mpc)  SigR(Msun/pc^2) err(weights) err(resampling) Nsources       <R>          e2\n");
         sign = -1.0;
      }
      for(i=0;i<para.nbins;i++){
         /* reminder: non-resampled value are stored from i=0 to nbins - 1 in result.w and result.GG */

         /* R and Rmean (weighted) */
         if(para.log){
            R     = meanR =exp(para.min+para.Delta*(double)i+para.Delta/2.0);
            meanR = exp(result.meanR[i]/result.w[i]);
         }else{
            R     = para.min+para.Delta*(double)i+para.Delta/2.0;
            meanR = result.meanR[i]/result.w[i];
         }

         /* don't divide by zero if there are too few objects in the bin */
         if(result.Nsources[i] > 0.0){
            fprintf(fileOut,"%12.7f %12.7f %12.7f %12.7f %15zd %12.7f %12.7f\n",
            R, sign*result.GG[i]/result.w[i],				\
            sqrt(1.0/result.w[i]), 	err_r[i],			\
            (long)result.Nsources[i],				\
            meanR,							\
            sign*result.e2[i]/result.w[i]);
         }else{
            fprintf(fileOut,"%12.7f %12.7f %12.7f %12.7f %15zd %12.7f %12.7f\n",
            R, 0.0,	0.0, 0.0, (long)result.Nsources[i], R, 0.0);
         }
      }
      fclose(fileOut);

      if(para.printSamples){

         /* TODO */
         /*
         sprintf(fileOutName, "%s.samples", para.fileOutName);
         fileOut = fopen(fileOutName, "w");


         fclose(fileOut);
         */
      }

      /*    covariance matrix */
      if(para.cov_mat){
         double *cov   = (double *)malloc(para.nbins*para.nbins*sizeof(double));
         sprintf(fileOutName, "%s.cov", para.fileOutName);
         fileOut = fopen(fileOutName, "w");
         for(i=0;i<para.nbins*para.nbins;i++) cov[i] = 0.0;
         for(i=0;i<para.nbins;i++){
            for(j=0;j<para.nbins;j++){
               for(l=0;l<para.nsamples;l++){
                  cov[para.nbins*i+j] += norm*(GG_mean[i] + result.GG[para.nbins*(l+1)+i]/result.w[para.nbins*(l+1)+i])*(GG_mean[j] + result.GG[para.nbins*(l+1)+j]/result.w[para.nbins*(l+1)+j]);
               }
               /*    write outfile */
               if(result.Nsources[i] > 3.0 && result.Nsources[j] > 3.0){
                  fprintf(fileOut,"%g ", cov[para.nbins*i+j]);
               }else{
                  fprintf(fileOut,"%f ", 0.0);
               }
            }
            fprintf(fileOut,"\n");
         }
         fclose(fileOut);
         free(cov);
      }
   }

   freeResult(para, result);

   return;
}
