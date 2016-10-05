#include "main.h"

/*
 *    main.c
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
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
 *    - when resampling, some objects fall outside the excluded zone
 *    - the mean <R> is R for auto_wp (why?)
 *    - errors for w(R) when using weights are zero
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
         /* number counts */
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

/*----------------------------------------------------------------*
*Main routines                                                   *
*----------------------------------------------------------------*/

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
   resample(&para, &random, dimStart, &mask, FIRSTCALL);

   /*    send data */
   comment(para, "sending data...");
   comData(para, &random, 0, dimStart, FIRSTCALL);
   comData(para, &data  , 0, dimStart, FIRSTCALL);

   /*    grow trees */
   comment(para, "building trees...");
   dataTree   = buildTree(&para, &data, &mask, dimStart, FIRSTCALL);     freePoint(para, data);
   comment(para, "done.\n");

   /*    Output tree to ascii file. format: RA DEC [w_0...w_nsamples] rank*/
   if(para.rank == 0 && para.printTree){
      sprintf(fileOutName, "%s.data.tree", para.fileOutName);
      printTree(para, fileOutName, dataTree, ROOT, 1, FIRSTCALL);
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
   resample(&para, &random, dimStart, &mask, FIRSTCALL);

   /*    send data */
   comment(para, "sending data...");
   comData(para, &random, 0, dimStart, FIRSTCALL);
   comData(para, &data  , 0, dimStart, FIRSTCALL);

   /*    grow trees */
   comment(para, "building trees...");
   randomTree = buildTree(&para, &random, &mask, dimStart, FIRSTCALL);   freePoint(para, random);
   dataTree   = buildTree(&para, &data, &mask, dimStart, FIRSTCALL);     freePoint(para, data);
   comment(para, "done.\n");

   /*    divide and conquer */
   long nodeSlaveRan  = splitTree(&para, &randomTree, ROOT, para.size, FIRSTCALL);
   long nodeSlaveData = splitTree(&para, &dataTree,   ROOT, para.size, FIRSTCALL);

   /*    output tree to ascii file. format: RA DEC [w_0...w_nsamples] rank*/
   if(para.rank == 0 && para.printTree){
      sprintf(fileOutName, "%s.data.tree", para.fileOutName);
      printTree(para, fileOutName, dataTree, ROOT, 1, FIRSTCALL);
      //exit(-1);
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
         comment(para, "RR(rp,pi)...       "); RR = Npairs3D(&para, &randomTree, ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
      }
      comment(para, "DR(rp,pi)...       "); DR = Npairs3D(&para, &dataTree,   ROOT, &randomTree, nodeSlaveRan,  FIRSTCALL);
      comment(para, "DD(rp,pi)...       "); DD = Npairs3D(&para, &dataTree,   ROOT, &dataTree,   nodeSlaveData, FIRSTCALL);
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
            fwrite(RR.NN,    sizeof(double), para.nbins*para.nbins*(para.nsamples+1), fileRR);
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
            RR.NN    = (double *)malloc(para.nbins*para.nbins*(para.nsamples+1)*sizeof(double));
            RR.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
            RR.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
            RR.meanR = (double *)malloc(para.nbins*sizeof(double));
            RR.NN_s  = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));

            fread(RR.NN,    sizeof(double), para.nbins*para.nbins*(para.nsamples+1), fileRR);
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
            for(j=0;j<para.nbins;j++) sum += DD.NN[i + para.nbins*j];
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
               DD_sum = 0;
               RR_sum = 0;
               j      = 0;
               while(R[j] < para.pi_max && j < para.nbins){ /* sum over pi integration */
                  DD_sum += DD.NN[i + para.nbins*j];
                  RR_sum += RR.NN[i + para.nbins*j];
                  j++;
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
               while(R[j] < para.pi_max && j < para.nbins){  /* sum over pi integration - ONLY informative */
                  RR_sum += RR.NN[i + para.nbins*j];
                  DR_sum += DR.NN[i + para.nbins*j];
                  DD_sum += DD.NN[i + para.nbins*j];
                  j++;
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
            for(j=0;j<para.nbins;j++){
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
   resample(&para, &random1, dimStart, &mask, FIRSTCALL);

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
   comment(para, "done.\n");

   /*    Ouput tree to ascii file. format: RA DEC [w_0...w_nsamples] rank */
   if(para.rank == 0 && para.printTree){

      sprintf(fileOutName, "%s.data1.tree", para.fileOutName);
      printTree(para, fileOutName, dataTree1, ROOT, 1, FIRSTCALL);

      sprintf(fileOutName, "%s.data2.tree", para.fileOutName);
      printTree(para, fileOutName, dataTree2, ROOT, 1, FIRSTCALL);
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
            comment(para, "R1R2(rp,pi)...       "); R1R2 = Npairs3D(&para, &randomTree1, ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         }

         // DEBUGGING - which one is correct for Peebles? D1R2 or D1R2?
         // comment(para, "D1R2(rp,pi)...       "); D1R2 = Npairs3D(&para, &dataTree1,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
         // comment(para, "D2R1(rp,pi)...       "); D2R1 = Npairs3D(&para, &dataTree2,   ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         comment(para, "D1R2(rp,pi)...       "); D1R2 = Npairs3D(&para, &dataTree1,   ROOT, &randomTree2, nodeSlaveRan2,  FIRSTCALL);
         comment(para, "D2R1(rp,pi)...       "); D2R1 = Npairs3D(&para, &dataTree2,   ROOT, &randomTree1, nodeSlaveRan1,  FIRSTCALL);
         comment(para, "D1D2(rp,pi)...       "); D1D2 = Npairs3D(&para, &dataTree1,   ROOT, &dataTree2,   nodeSlaveData2, FIRSTCALL);
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
            fwrite(R1R2.NN,    sizeof(double), para.nbins*para.nbins*(para.nsamples+1), fileRR);
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
            R1R2.NN    = (double *)malloc(para.nbins*para.nbins*(para.nsamples+1)*sizeof(double));
            R1R2.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
            R1R2.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
            R1R2.meanR = (double *)malloc(para.nbins*sizeof(double));
            R1R2.NN_s  = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));

            fread(R1R2.NN,    sizeof(double), para.nbins*para.nbins*(para.nsamples+1), fileRR);
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
            for(j=0;j<para.nbins;j++) sum += D1D2.NN[i + para.nbins*j];
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
               while(R[j] < para.pi_max && j < para.nbins){ /* sum over pi integration */
                  D1D2_sum += D1D2.NN[i + para.nbins*j];
                  R1R2_sum += R1R2.NN[i + para.nbins*j];
                  j++;
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
               while(R[j] < para.pi_max && j < para.nbins){  /* sum over pi integration */
                  R1R2_sum += R1R2.NN[i + para.nbins*j];
                  D1R2_sum += D1R2.NN[i + para.nbins*j];
                  D2R1_sum += D2R1.NN[i + para.nbins*j];
                  D1D2_sum += D1D2.NN[i + para.nbins*j];
                  j++;
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
            for(j=0;j<para.nbins;j++){
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
   resample(&para, &source, dimStart, &mask, FIRSTCALL);

   /*    send data. Source catalogue is partitioned among cpus */
   comment(para, "sending data...");

   comData(para, &source, para.size, dimStart, FIRSTCALL);
   comData(para, &lens  , 0, dimStart, FIRSTCALL);

   /*    grow trees */
   comment(para, "building trees...");
   sourceTree = buildTree(&para, &source, &mask, dimStart, FIRSTCALL);   freePoint(para, source);
   lensTree   = buildTree(&para, &lens,   &mask, dimStart, FIRSTCALL);   freePoint(para, lens);

   /*    TODO */
   if(para.rank == 0 && para.printTree){
      if(para.verbose) fprintf(stderr,"\n%s: **WARNING** printTree not supported for gglens.\n", MYNAME);
   }

   comment(para, "done.\n");
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

double wTheta(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R2, Result D2R1, int i, int l){
   /* i is the bin index. l is sample 0 to 256 (0: no resampling, 1 -> 256: bootstrap or jackknife samples) ) */

   /*    initialization */
   // double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]-1))/(D1D2.N1[l]*(D1D2.N2[l]-1));
   // double Norm2 = (R1R2.N2[l]-1)/D1D2.N1[l];
   // double Norm3 = (R1R2.N1[l]-1)/D1D2.N2[l];
   double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]))/(D1D2.N1[l]*(D1D2.N2[l]));
   double Norm2 = (R1R2.N1[l])/D1D2.N1[l];
   double Norm3 = (R1R2.N2[l])/D1D2.N2[l];

   double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/((R1R2.N2[l]-1)*(D1D2.N2[l]-1));
   double Norm5 = R1R2.N1[l]/(D1D2.N2[l]-1.0);


   double result = 0.0;

   if(D1D2.NN[para.nbins*l+i] > 0
      && D1R2.NN[para.nbins*l+i] > 0
      && D2R1.NN[para.nbins*l+i] > 0
      && R1R2.NN[para.nbins*l+i] > 0){
         switch(estimator){
            case LS:  /*   Landy and Szalay */
            result  =  Norm1*D1D2.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i];
            result += -Norm2*D1R2.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i];
            result += -Norm3*D2R1.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i] + 1.0;
            break;
            case NAT: /*   Natural */
            result = Norm1*D1D2.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i] - 1.0;
            break;
            case HAM: /*   Hamilton */
            result = Norm4*D1D2.NN[para.nbins*l+i]*R1R2.NN[para.nbins*l+i]/(D1R2.NN[para.nbins*l+i]*D2R1.NN[para.nbins*l+i]) - 1.0;
            break;
            case PEEBLES: /*  Peebles */
            result = Norm5*D1D2.NN[para.nbins*l+i]/D1R2.NN[para.nbins*l+i] - 1.0;
            break;
         }
      }

      return result;
   }

   double wp(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R2, Result D2R1, int i, int j, int l){
      /* i,j are the bin indexes. i : pi, j : rp.
      if i = -1, integrates over i (pi).
      l is sample 0 to 256 (0: no resampling, 1 -> 256: bootstrap or jackknife samples) ) */
      double R;

      /*    initialization */
      // double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]-1))/(D1D2.N1[l]*(D1D2.N2[l]-1));
      // double Norm2 = (R1R2.N2[l]-1)/D1D2.N1[l];
      // double Norm3 = (R1R2.N1[l]-1)/D1D2.N2[l];
      double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]))/(D1D2.N1[l]*(D1D2.N2[l]));
      double Norm2 = (R1R2.N1[l])/D1D2.N1[l];
      double Norm3 = (R1R2.N2[l])/D1D2.N2[l];


      double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/((R1R2.N2[l]-1)*(D1D2.N2[l]-1));
      double Norm5 = R1R2.N1[l]/(D1D2.N2[l]-1.0);

      double result = 0.0, sum = 0.0;

      if(i > -1){
         // if( D1D2.NN[j + para.nbins*(i + para.nbins*l)] > 0
         //&& D1R2.NN[j + para.nbins*(i + para.nbins*l)] > 0
         //&& D2R1.NN[j + para.nbins*(i + para.nbins*l)] > 0
         //&& R1R2.NN[j + para.nbins*(i + para.nbins*l)] > 0){

         if(R1R2.NN[j + para.nbins*(i + para.nbins*l)] > 0){
            switch(estimator){
               case LS:  /*   Landy and Szalay */
               result  =  Norm1*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)];
               result += -Norm2*D1R2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)];
               result += -Norm3*D2R1.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)] + 1.0;
               break;
               case NAT: /*    Natural */
               result = Norm1*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)] - 1.0;
               break;
               case HAM: /*   Hamilton */
               result = Norm4*D1D2.NN[j + para.nbins*(i + para.nbins*l)]*R1R2.NN[j + para.nbins*(i + para.nbins*l)]/(D1R2.NN[j + para.nbins*(i + para.nbins*l)]*D2R1.NN[j + para.nbins*(i + para.nbins*l)]) - 1.0;
               break;
               case PEEBLES: /*  Peebles */
               result = Norm5*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/D1R2.NN[j + para.nbins*(i + para.nbins*l)] - 1.0;
               break;
            }
         }
         return result;
      }else{
         /*    reset i and integrate over pi (i)*/
         // i = 0;
         // R = 0.0;
         result = 0.0;

         // while(R < para.pi_max && i < para.nbins){
         for(i=0; i<para.nbins; i++){

            sum = 0.0;
            //     if( D1D2.NN[j + para.nbins*(i + para.nbins*l)] > 0
            //	  && D1R2.NN[j + para.nbins*(i + para.nbins*l)] > 0
            //	  && D2R1.NN[j + para.nbins*(i + para.nbins*l)] > 0
            //	  && R1R2.NN[j + para.nbins*(i + para.nbins*l)] > 0){

            if( R1R2.NN[j + para.nbins*(i + para.nbins*l)] > 0.0){
               switch(estimator){
                  case LS:  /* Landy and Szalay */
                  sum  =  Norm1*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)];
                  sum += -Norm2*D1R2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)];
                  sum += -Norm3*D2R1.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)] + 1.0;
                  break;
                  case NAT: /* Natural */
                  sum = Norm1*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins*l)] - 1.0;
                  break;
                  case HAM: /* Hamilton */
                  sum = Norm4*D1D2.NN[j + para.nbins*(i + para.nbins*l)]*R1R2.NN[j + para.nbins*(i + para.nbins*l)]/(D1R2.NN[j + para.nbins*(i + para.nbins*l)]*D2R1.NN[j + para.nbins*(i + para.nbins*l)]) - 1.0;
                  break;
                  case PEEBLES: /* Peebles */
                  sum = Norm5*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/D1R2.NN[j + para.nbins*(i + para.nbins*l)] - 1.0;
                  break;
               }
            }

            result += 2.0*sum*para.Delta_pi;

            // i++;
            /* R */
            // R = 0.0+para.Delta_pi*(double)i;

         }
         return result;
      }
   }


   double xis(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R2, Result D2R1, int i, int l){
      /*    i is the bin index. l is sample 0 to 256 (0: no resampling, 1 -> 256: bootstrap or jackknife samples) ) */

      /*    initialization */
      // double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]-1))/(D1D2.N1[l]*(D1D2.N2[l]-1));
      // double Norm2 = (R1R2.N2[l]-1)/D1D2.N1[l];
      // double Norm3 = (R1R2.N1[l]-1)/D1D2.N2[l];
      double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]))/(D1D2.N1[l]*(D1D2.N2[l]));
      double Norm2 = (R1R2.N1[l])/D1D2.N1[l];
      double Norm3 = (R1R2.N2[l])/D1D2.N2[l];



      double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/((R1R2.N2[l]-1)*(D1D2.N2[l]-1));
      double Norm5 = R1R2.N1[l]/(D1D2.N2[l]-1.0);

      double result = 0.0;

      if(D1D2.NN_s[para.nbins*l+i] > 0
         && D1R2.NN_s[para.nbins*l+i] > 0
         && D2R1.NN_s[para.nbins*l+i] > 0
         && R1R2.NN_s[para.nbins*l+i] > 0){
            switch(estimator){
               case LS:  /*   Landy and Szalay */
               result  =  Norm1*D1D2.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i];
               result += -Norm2*D1R2.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i];
               result += -Norm3*D2R1.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i] + 1.0;
               break;
               case NAT: /*   Natural */
               result = Norm1*D1D2.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i] - 1.0;
               break;
               case HAM: /*   Hamilton */
               result = Norm4*D1D2.NN_s[para.nbins*l+i]*R1R2.NN_s[para.nbins*l+i]/(D1R2.NN_s[para.nbins*l+i]*D2R1.NN_s[para.nbins*l+i]) - 1.0;
               break;
               case PEEBLES: /*  Peebles */
               result = Norm5*D1D2.NN_s[para.nbins*l+i]/D1R2.NN_s[para.nbins*l+i] - 1.0;
               break;
            }
         }

         return result;
      }

      #define leaf(tree,node) ((tree)->left[(node)] == -1)
      #define node(tree,node) ((tree)->left[(node)] > -1)

      Result Nobjects(const Config *para, const Tree *tree1, const long i, int firstCall){
         /*    returns the number of objects in tree1.
         *     Uses random objects to cut the
         *     area into constant parts.
         */

         long k = 0, l;
         static Result result;
         static long count, total;
         double w, NN;

         if(firstCall){

            count = 0;
            total = tree1->N[i];

            /*    number counts */
            result.NN = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));
            for(k = 0; k < para->nbins*(para->nsamples+1); k++) result.NN[k] = 0;

            /*    weighted center of the bin */
            result.meanR = (double *)malloc(para->nbins*sizeof(double));
            for(k = 0; k < para->nbins; k++) result.meanR[k]    = 0.0;

            /*    total number of objects per sample */
            result.N1 = (double *)malloc((para->nsamples+1)*sizeof(double));
            result.N2 = (double *)malloc((para->nsamples+1)*sizeof(double));
            result.w  = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));

            result.N1[0] = tree1->N[ROOT];
            for(l=0;l<para->nsamples;l++){
               result.N1[l+1] = tree1->Ntot[l];
            }
         }

         /*    compute the number of objects recursively */
         if(node(tree1, i)){

            Nobjects(para, tree1, tree1->left[i],  0);
            Nobjects(para, tree1, tree1->right[i], 0);

         }else{

            NN = tree1->point.x[NDIM*i+2];
            if(para->log) NN = log(NN);

            k = floor((NN - para->min)/para->Delta);
            if(0 <= k && k < para->nbins){
               if(para->weighted) w = tree1->point.w[i];
               else w = 1.0;

               result.NN[k] +=  1.0*w;
               for(l=0;l<para->nsamples;l++){

               /* DEBUGGING for SUBSAMPLE
               double weight =  tree1->w[para->nsamples*i + l];
               if(weight == 0){
                  weight = 1;
               }else{
                  weight = 0;
               }
               result.NN[para->nbins*(l+1) + k] += weight;
               */
            result.NN[para->nbins*(l+1) + k] += tree1->w[para->nsamples*i + l]*w;
         }
         result.meanR[k] += NN*w;
      }
      count += tree1->N[i];
      printCount(count, total, 1000, para->verbose);

   }

   if(firstCall && para->verbose) fprintf(stderr, "\b\b\b\b\b\b\b%6.2f%%\n", 100.0);

   return result;
}


#define z_tree1         ((tree1)->point.x[NDIM*(i)+2])
#define distComo_tree1  ((tree1)->distComo[(i)])
#define z_tree2         ((tree2)->point.x[NDIM*(i)+2])
#define distComo_tree2  ((tree2)->distComo[(i)])


Result Npairs(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall){
   /*    returns the number of pairs if tree1 != tree2, and twice
    *    the number of pairs if tree1 == tree2, so that the treatment of pairs
    *    in w(theta) estimators is identical (numerical trick).
    */

   long k = 0, l;
   static Result result;
   static long count, total;
   double NN, d, deltaTheta;

   if(firstCall){
      count = 0;

      result.NN = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));
      result.w  = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));
      for(k = 0; k < para->nbins*(para->nsamples+1); k++){
         result.NN[k] = 0;
         result.w[k]  = 0.0;
      }
      result.meanR    = (double *)malloc(para->nbins*sizeof(double));
      for(k = 0; k < para->nbins; k++){
         result.meanR[k]    = 0.0;
      }

      total = tree1->N[i]*tree2->N[j];

      /* number of points (used in wTheta(...) for the normalization) */
      result.N1 = (double *)malloc((para->nsamples+1)*sizeof(double));
      result.N2 = (double *)malloc((para->nsamples+1)*sizeof(double));
      result.N1[0] = tree1->N[ROOT];
      result.N2[0] = tree2->N[ROOT];

      for(l=0;l<para->nsamples;l++){
         result.N1[l+1] = tree1->Ntot[l];
         result.N2[l+1] = tree2->Ntot[l];

      }
   }

   if((tree1 == tree2 && i > j) || j < 0) return result;

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

      if(para->corr != AUTO_3D && para->corr != CROSS_3D){
         switch(para->proj){
            case COMO:
               d  = distComo_tree1*deltaTheta*PI/180.0;              /*    Transverse distance in comoving coordinates (Mpc) */
               break;
            case PHYS:
               d  = distComo_tree1*deltaTheta*PI/180.0/(1+z_tree1);  /*    Transverse distance in physical coordinates (Mpc) */
               break;
            case THETA:
               d = deltaTheta;
               break;
         }
      }else{
         d = deltaTheta;
      }


      /*    Note: tree->N[i] is the weighted sum of objects so there's no
       *    need to keep track of the sum of the weights, which is simply NN.
       */
      if(para->log) d = log(d);
      k = floor((d - para->min)/para->Delta);
      if(0 <= k && k < para->nbins){
         NN = tree1->N[i]*tree2->N[j];
         result.NN[k] += NN;
         for(l=0;l<para->nsamples;l++){
            result.NN[para->nbins*(l+1) + k] += NN*tree1->w[para->nsamples*i + l]*tree2->w[para->nsamples*j + l];
         }
         result.meanR[k] += d*NN;
      }
      count += tree1->N[i]*tree2->N[j];
      printCount(count,total,10000,para->verbose);
   }

   if(firstCall && para->verbose) fprintf(stderr, "\b\b\b\b\b\b\b%6.2f%%\n", 100.0);

   if(firstCall && tree1 == tree2){
      for(k = 0; k < para->nbins*(para->nsamples+1); k++){
         result.NN[k] *= 2;
      }
      for(k = 0; k < para->nbins; k++){
         result.meanR[k] *= 2.0;
      }
   }

   return result;
}

Result Npairs3D(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall){
   /*    returns the number of pairs if tree1 != tree2, and twice
    *    the number of pairs if tree1 == tree2, so that the treatment of pairs
    *    in w(theta) estimators is identical (numerical trick).
    */

   long k = 0, m = 0, p = 0, l;
   static Result result;
   static long count, total;
   double NN, deltaTheta, rp, pi, s, d;

   if(firstCall){
      count = 0;
      total = tree1->N[i]*tree2->N[j];
      result.NN = (double *)malloc(para->nbins*para->nbins*(para->nsamples+1)*sizeof(double));
      for(k = 0; k < para->nbins*para->nbins*(para->nsamples+1); k++) result.NN[k] = 0.0;
      /*    for xi(s) */
      result.NN_s = (double *)malloc(para->nbins*(para->nsamples+1)*sizeof(double));
      for(p = 0; p < para->nbins*(para->nsamples+1); p++) result.NN_s[p] = 0.0;

      result.meanR = (double *)malloc(para->nbins*sizeof(double));
      for(k = 0; k < para->nbins; k++){
         result.meanR[k]    = 0.0;
      }
      /*    number of points (used in wTheta(...) for the normalization) */
      result.N1 = (double *)malloc((para->nsamples+1)*sizeof(double));
      result.N2 = (double *)malloc((para->nsamples+1)*sizeof(double));
      result.N1[0] = tree1->N[ROOT];
      result.N2[0] = tree2->N[ROOT];
      for(l=0;l<para->nsamples;l++){
         result.N1[l+1] = tree1->Ntot[l];
         result.N2[l+1] = tree2->Ntot[l];
      }
   }

   if((tree1 == tree2 && i > j) || j < 0) return result;

   deltaTheta = para->distAng(tree1, &i, tree2, &j);

   // DEBUGGING for wp(r) testing
   //deltaTheta = distAngCart(tree1, &i, tree2, &j);

   d = deltaTheta;

   if(node(tree1, i) && tree1->r[i]/d > para->OA){
      if(node(tree2, j) && tree2->r[j]/d > para->OA){
         Npairs3D(para, tree1, tree1->left[i],  tree2, tree2->left[j],  0);
         Npairs3D(para, tree1, tree1->right[i], tree2, tree2->left[j],  0);
         Npairs3D(para, tree1, tree1->left[i],  tree2, tree2->right[j], 0);
         Npairs3D(para, tree1, tree1->right[i], tree2, tree2->right[j], 0);
      }else{
         Npairs3D(para, tree1, tree1->left[i], tree2,  j,  0);
         Npairs3D(para, tree1, tree1->right[i], tree2, j,  0);
      }
   }else if(node(tree2,j) && tree2->r[j]/d > para->OA){
      Npairs3D(para, tree1, i, tree2, tree2->left[j],   0);
      Npairs3D(para, tree1, i, tree2,  tree2->right[j], 0);
   }else{

      pi = ABS(tree1->distComo[i] - tree2->distComo[j]);
      rp = (tree1->distComo[i]+tree2->distComo[j])/2.0*deltaTheta*PI/180.0;

      // DEBUGGING for wp(r) testing
      //double delx = (tree1->point.x[NDIM*(i)+0] - tree2->point.x[NDIM*(j)+0]);
      //double dely = (tree1->point.x[NDIM*(i)+1] - tree2->point.x[NDIM*(j)+1]);
      //double delz = (tree1->point.x[NDIM*(i)+2] - tree2->point.x[NDIM*(j)+2]);
      //rp = sqrt((delx*delx)+(dely*dely));
      //pi = ABS(delz);

      s  =  sqrt(pi*pi + rp*rp);

      if(para->proj == PHYS){      /* Distances in physical coordinates (Mpc) */
         pi /= 1.0+z_tree1;
         rp /= 1.0+z_tree1;
         s  /= 1.0+z_tree1;
      }

      if(para->log){
         rp = log(rp);
         s  = log(s);
      }

      /* pi in linear scale */
      k = floor((pi - 0.0)/para->Delta_pi);
      m = floor((rp - para->min)/para->Delta);

      if(0 <= k && k < para->nbins && 0 <= m && m < para->nbins){
         NN = tree1->N[i]*tree2->N[j];
         result.NN[m + para->nbins*k] += NN;
         for(l=0;l<para->nsamples;l++){
            result.NN[m + para->nbins*(k + (para->nbins*(l+1)))] += NN*tree1->w[para->nsamples*i + l]*tree2->w[para->nsamples*j + l];
         }
         result.meanR[m] += rp*NN;
      }

      p = floor((s  - para->min)/para->Delta);

      if(0 <= p && p < para->nbins){
         NN = tree1->N[i]*tree2->N[j];
         result.NN_s[p] += NN;
         for(l=0;l<para->nsamples;l++){
            result.NN_s[p + (para->nbins*(l+1))] += NN*tree1->w[para->nsamples*i + l]*tree2->w[para->nsamples*j + l];
         }
      }

      count += tree1->N[i]*tree2->N[j];
      printCount(count, total, 10000, para->verbose);
   }

   if(firstCall && para->verbose) fprintf(stderr, "\b\b\b\b\b\b\b%6.2f%%\n",100.0);

   if(firstCall && tree1 == tree2){
      for(k = 0; k < para->nbins*para->nbins*(para->nsamples+1); k++){
         result.NN[k] *= 2;
      }
      for(p = 0; p < para->nbins*(para->nsamples+1); p++){
         result.NN_s[p] *= 2;
      }
      for(k = 0; k < para->nbins; k++){
         result.meanR[k] *= 2;
      }

   }

   return result;
}

#undef z_tree2
#undef distComo_tree2
#undef z_tree1
#undef distComo_tree1

Result gg(const Config *para, const Tree *lens, const long i, const Tree *source, const long j, int firstCall){
   /*    Computes the galaxy-galaxy two-point correlation function. */

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

   /*    correlate if (leaf OR size/d < OA), otherwise go down the tree. Redshift is point.x[NDIM*i + 2] */
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
   /*    See Leauthaud et al. (2010),  2010ApJ...709...97L */

   double dA, dR, invScaleFac;
   long l, k, zero = 0;

   /*    very quick tests? But might be time consuming, too. The key
    *    is to find the best balance...
    *    e.g. if(z_source < z_lens) return;
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

   /*    dA = phy_dis/angular_size_in_radians = tran_como_dis/(1+z)*/
   dA  = distComo_lens/(1.0 + z_lens);               /*  Angular diameter distance in physical coordinates */
   dR  = dA*deltaTheta*PI/180.0;                     /*  Transverse distance in phys coordinates (Mpc)     */
   dR *= invScaleFac;                                /*  If coordinate in comoving system, divide by a     */

   if(para->log) dR = log(dR);
   k = floor((dR - para->min)/para->Delta);
   if(0 <= k && k < para->nbins && z_source > z_lens + zerr_lens + zerr_source && z_source > zerr_lens + para->deltaz){

      double DOS        = distComo_source;                  /*    this is (tranverse) comoving distance   */
      double DOL        = distComo_lens;                    /*    but (1+z_source) factors cancel out     */
      double DLS        = distComo_source - distComo_lens;  /*    Approx. Omega_k = 0                      */
      double SigCritInv = DOL*DLS/DOS/1000.0/(1.0 + z_lens);/*    1/SigCrit in Gpc (1+z) to convert to dA  */
      SigCritInv       /= 1.663e3;                          /*    see Narayan & Bartelmann pge 10          */
      SigCritInv       *= invScaleFac*invScaleFac;          /*    If coordinates in comoving system        */

      /*    lensing weight */
      double WW, GG, w   = SigCritInv*SigCritInv*w_source;

      if(para->calib){ /*  Compute calibration factor (1+m or c) */

         GG = e1_source*w/SigCritInv;
         WW = w/SigCritInv;
         //WW = w;

      }else{           /*  Compute gg lensing signal */

         /*    Point A */
         Point A = createPoint(*para, 1);
         A.x[0]    = RA_source;
         A.x[1]    = DEC_lens;
         double AS = distAngPointSpher(para, &A, &zero, &source->point, &j);
         double AL = distAngPointSpher(para, &A, &zero, &lens->point,   &i);
         freePoint(*para, A);

         /*       to get correct sign for phi_gg */
         if(RA_source  > RA_lens)  AL = -AL;
         if(DEC_source < DEC_lens) AS = -AS;

         double phi_gg     = atan2(AS, AL);
         double cos2phi_gg = cos(2.0*phi_gg);
         double sin2phi_gg = sin(2.0*phi_gg);
         double e1         =  e1_source*cos2phi_gg + e2_source*sin2phi_gg;
         double e2         = -e1_source*sin2phi_gg + e2_source*cos2phi_gg;

         GG = e1*w/SigCritInv;
         WW = w;

         result.e2[k] += e2*w/SigCritInv;

      }

      /*    signal */
      result.GG[k] += GG;
      result.w[k]  += WW;

      /*    bootstraps */
      for(l=0;l<para->nsamples;l++){
         /*    ATTENTION: source->w and lens->w below are resampling weights (boostrap or jackknife) whereas
          *    w above (namely result.w) is the lensing weight (i.e. shape measurement error
          *    if inverse variance estimate)
          */
         result.GG[para->nbins*(l+1) + k] += GG*lens->w[para->nsamples*i + l]*source->w[para->nsamples*j + l];
         result.w[para->nbins*(l+1) + k]  += WW*lens->w[para->nsamples*i + l]*source->w[para->nsamples*j + l];
      }

      /*    keep track of info per bin */
      result.Nsources[k] += 1.0;
      result.meanR[k]    += dR*WW;
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

#undef leaf
#undef node

void freeResult(const Config para, Result result){

   switch(para.corr){
      case AUTO: case CROSS: case AUTO_3D: case CROSS_3D: case AUTO_WP: case CROSS_WP: case NUMBER:
      free(result.N1);
      free(result.N2);
      free(result.NN);
      free(result.meanR);
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



/*
 *    MPI routines
 */

#define BASE 200
void comData(const Config para, Point *data, long Ncpu, int dim, int firstCall){
   /*    Partitions and sends the data recursively if para.rank == master.
    *    Otherwise receives the data. Then returns the final data
    */

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

         /*    next splitting coordinate */
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
            MPI_Send(data->w   , data->N, MPI_DOUBLE, rank, BASE+6, MPI_COMM_WORLD);
            if(para.corr == GGLENS){
               MPI_Send(data->zerr, data->N, MPI_DOUBLE, rank, BASE+3, MPI_COMM_WORLD);
               MPI_Send(data->e1  , data->N, MPI_DOUBLE, rank, BASE+4, MPI_COMM_WORLD);
               MPI_Send(data->e2  , data->N, MPI_DOUBLE, rank, BASE+5, MPI_COMM_WORLD);
            }

         }else{
            /*    master can't send data to itself, so we're using a little trick here.
             *    It first recovers the pointers from the first call ("Parent"), then
             *    shifts the data back to the start, and finally reallocates the array
             *    to free the unsused memory (i.e. data sent to other cpus).
             *    Note: memmove() is safer than memcpy().
             */

            dataParent->x = (double *)memmove(dataParent->x, data->x, data->N*NDIM*sizeof(double));
            dataParent->w    = (double *)memmove(dataParent->w,    data->w ,   data->N*sizeof(double));
            if(para.corr == GGLENS){
               dataParent->zerr = (double *)memmove(dataParent->zerr, data->zerr, data->N*sizeof(double));
               dataParent->e1   = (double *)memmove(dataParent->e1,   data->e1,   data->N*sizeof(double));
               dataParent->e2   = (double *)memmove(dataParent->e2,   data->e2,   data->N*sizeof(double));
            }
            dataParent->x = (double *)realloc(dataParent->x, data->N*NDIM*sizeof(double));
            dataParent->w    = (double *)realloc(dataParent->w,    data->N*sizeof(double));
            if(para.corr == GGLENS){
               dataParent->zerr = (double *)realloc(dataParent->zerr, data->N*sizeof(double));
               dataParent->e1   = (double *)realloc(dataParent->e1,   data->N*sizeof(double));
               dataParent->e2   = (double *)realloc(dataParent->e2,   data->N*sizeof(double));
            }

            /*    copy the new number of points and the dim
             *    along which they are sorted
             */
            dataParent->N   = data->N;
            dataParent->dim = dim;
         }

         rank--;
         /*    call recursively to transmit the entire data set to all cpus without further splitting */
         if(Ncpu == 0 && rank >= 0) comData(para, data, 0, dim, 0);
      }

   }else{
      MPI_Status status;

      MPI_Recv(&data->N, 1, MPI_LONG, MASTER, BASE+0, MPI_COMM_WORLD, &status);
      MPI_Recv(&dim, 1, MPI_INT, MASTER, BASE+1, MPI_COMM_WORLD, &status);

      data->x = (double *)malloc(data->N*NDIM*sizeof(double));
      MPI_Recv(data->x, data->N*NDIM, MPI_DOUBLE, MASTER, BASE+2, MPI_COMM_WORLD, &status);

      data->w    = (double *)malloc(data->N*sizeof(double));
      MPI_Recv(data->w   , data->N, MPI_DOUBLE, MASTER, BASE+6, MPI_COMM_WORLD, &status);

      if(para.corr == GGLENS){
         data->zerr = (double *)malloc(data->N*sizeof(double));
         data->e1   = (double *)malloc(data->N*sizeof(double));
         data->e2   = (double *)malloc(data->N*sizeof(double));
         MPI_Recv(data->zerr, data->N, MPI_DOUBLE, MASTER, BASE+3, MPI_COMM_WORLD, &status);
         MPI_Recv(data->e1  , data->N, MPI_DOUBLE, MASTER, BASE+4, MPI_COMM_WORLD, &status);
         MPI_Recv(data->e2  , data->N, MPI_DOUBLE, MASTER, BASE+5, MPI_COMM_WORLD, &status);
      }
   }

   return;
}
#undef BASE

#define BASE 300
void comResult(const Config para, Result result, long Ncpu, int split){
   /*    Sends results back to master if slave. Sums up everything if master. If
    *    split is larger than 0, it means the data has been split before building
    *    the tree, hence it requires to sum up the total number of objects as well,
    *    but only N2, because N1 is never partitioned.
    */
   long i, rank;

   if(para.rank != MASTER){
      switch(para.corr){
         case AUTO: case CROSS: case AUTO_3D: case CROSS_3D: case NUMBER:
         MPI_Send(result.NN, para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+0, MPI_COMM_WORLD);
         MPI_Send(result.N1, para.nsamples+1, MPI_DOUBLE, MASTER, BASE+1, MPI_COMM_WORLD);
         MPI_Send(result.N2, para.nsamples+1, MPI_DOUBLE, MASTER, BASE+2, MPI_COMM_WORLD);
         MPI_Send(result.meanR,   para.nbins, MPI_DOUBLE, MASTER, BASE+3, MPI_COMM_WORLD);
         break;
         case GGLENS:
         MPI_Send(result.GG, para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+0, MPI_COMM_WORLD);
         MPI_Send(result.w,  para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+1, MPI_COMM_WORLD);
         MPI_Send(result.Nsources, para.nbins, MPI_DOUBLE, MASTER, BASE+2, MPI_COMM_WORLD);
         MPI_Send(result.meanR,    para.nbins, MPI_DOUBLE, MASTER, BASE+3, MPI_COMM_WORLD);
         MPI_Send(result.e2,       para.nbins, MPI_DOUBLE, MASTER, BASE+4, MPI_COMM_WORLD);
         break;
         case AUTO_WP: case CROSS_WP:
         MPI_Send(result.NN,    para.nbins*para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+0, MPI_COMM_WORLD);
         MPI_Send(result.N1,    para.nsamples+1, MPI_DOUBLE, MASTER, BASE+1, MPI_COMM_WORLD);
         MPI_Send(result.N2,    para.nsamples+1, MPI_DOUBLE, MASTER, BASE+2, MPI_COMM_WORLD);
         MPI_Send(result.meanR, para.nbins, MPI_DOUBLE, MASTER, BASE+3, MPI_COMM_WORLD);
         MPI_Send(result.NN_s,  para.nbins*(para.nsamples+1), MPI_DOUBLE, MASTER, BASE+4, MPI_COMM_WORLD);
         break;
      }
   }else{

      Result slave;
      MPI_Status status;

      switch(para.corr){
         case AUTO: case CROSS: case AUTO_3D: case CROSS_3D:  case NUMBER:
         slave.NN    = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(long));
         slave.N1    = (double *)malloc((para.nsamples+1)*sizeof(double));
         slave.N2    = (double *)malloc((para.nsamples+1)*sizeof(double));
         slave.meanR = (double *)malloc(para.nbins*sizeof(double));
         for(rank=1;rank<Ncpu;rank++){
            MPI_Recv(slave.NN,    para.nbins*(para.nsamples+1), MPI_DOUBLE, rank, BASE+0, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.N1,    para.nsamples+1, MPI_DOUBLE, rank, BASE+1, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.N2,    para.nsamples+1, MPI_DOUBLE, rank, BASE+2, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.meanR, para.nbins, MPI_DOUBLE, rank, BASE+3, MPI_COMM_WORLD, &status);
            for(i=0;i<para.nbins*(para.nsamples+1);i++) result.NN[i] += slave.NN[i];
            for(i=0;i<para.nbins;i++){
               result.meanR[i]    += slave.meanR[i];
            }
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
         case AUTO_WP: case CROSS_WP:
         slave.NN = (double *)malloc(para.nbins*para.nbins*(para.nsamples+1)*sizeof(double));
         slave.N1 = (double *)malloc((para.nsamples+1)*sizeof(double));
         slave.N2 = (double *)malloc((para.nsamples+1)*sizeof(double));
         slave.meanR = (double *)malloc(para.nbins*para.nbins*sizeof(double));
         slave.NN_s = (double *)malloc(para.nbins*(para.nsamples+1)*sizeof(double));
         for(rank=1;rank<Ncpu;rank++){
            MPI_Recv(slave.NN,    para.nbins*para.nbins*(para.nsamples+1), MPI_DOUBLE, rank, BASE+0, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.N1,    para.nsamples+1, MPI_DOUBLE, rank, BASE+1, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.N2,    para.nsamples+1, MPI_DOUBLE, rank, BASE+2, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.meanR, para.nbins, MPI_DOUBLE, rank, BASE+3, MPI_COMM_WORLD, &status);
            MPI_Recv(slave.NN_s,  para.nbins*(para.nsamples+1), MPI_DOUBLE, rank, BASE+4, MPI_COMM_WORLD, &status);
            for(i=0;i<para.nbins*para.nbins*(para.nsamples+1);i++) result.NN[i] += slave.NN[i];
            for(i=0;i<para.nbins*(para.nsamples+1);i++) result.NN_s[i] += slave.NN_s[i];
            for(i=0;i<para.nbins;i++){
               result.meanR[i] += slave.meanR[i];
            }
            if(split){
               /* only N2 may be partitioned */
               for(i=0;i<para.nsamples+1;i++) result.N2[i] += slave.N2[i];
            }
         }
         break;
         freeResult(para, slave);
      }
   }
}
#undef BASE
