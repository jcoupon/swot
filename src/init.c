#include "init.h"

/*
 *    init.c
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */



/*
 * 	Configuration
 */

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
   strcpy(MYNAME, "swot");

   /* default parameters below */
   strcpy(para->fileInName1, "data1.cat");
   strcpy(para->fileInName2, "data2.cat");
   strcpy(para->fileRanName1, "ran1.cat");
   strcpy(para->fileRanName2, "ran2.cat");
   strcpy(para->RRInFileName, "");
   strcpy(para->RROutFileName, "");
   /* colum ids, default:
   ids           1  2   3   4   5  6  7
   [lensing]     RA DEC z  deltaz e1 e2 weight
   [w(theta)]    RA DEC
   [wp(rp)]      RA DEC z
   [xi(r)]       X  Y   Z (Mpc)
   [number]      RA DEC X */


	//para->data1Id = NULL;
	//para->data2Id = NULL;
	//para->ran1Id = NULL;
	//para->ran2Id = NULL;

   for(i=0;i<NIDSMAX;i++){
		para->data1Id[i] = malloc((72+1) * sizeof(char));
		para->data2Id[i] = malloc((72+1) * sizeof(char));
		para->ran1Id[i] = malloc((72+1) * sizeof(char));
		para->ran2Id[i] = malloc((72+1) * sizeof(char));
      sprintf(para->data1Id[i], "%d", i+1);
      sprintf(para->data2Id[i], "%d", i+1);
      sprintf(para->ran1Id[i], "%d", i+1);
      sprintf(para->ran2Id[i], "%d", i+1);
      //para->data1Id[i] = i+1;
      //para->data2Id[i] = i+1;
      //para->ran1Id[i]  = i+1;
      //para->ran2Id[i]  = i+1;
   }
   NDIM              = 2; /* number of dimensions, 2 or 3 */
   para->distAng     = &distAngSpher;
   para->proj        = THETA;
   para->corr        = AUTO;
   para->estimator   = LS;
   para->min         = 0.0001;
   para->max         = 1.0;
   para->nbins       = 20;
   para->nbins_pi    = 100;
   para->log         = 1;
   para->OA          = 0.05;
   para->nsub        = 32;
   para->nsamples    = para->nsub;
   para->err         = JACKKNIFE;
   para->cov_mat     = 0;
   para->xi          = 0;
   para->deltaz      = 0.03;
   para->pi_max      = 40.0;
   para->weighted    = 0;
   strcpy(para->fileOutName,   "corr.out");
   para->calib       = 0;
   para->resample2D  = 0;
   para->printTree   = 0;
   para->printTreeAndExit = 0;
   para->seed         = -1;
   para->printSamples = 0;
   para->fits         = -1;

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
       (Super W Of Theta) MPI version 1.3\n\n\
Program to compute two-point correlation functions.\n\
Usage: %s -c configFile [options]: run the program\n\
       %s -d: display a default configuration file\n\
Important: for RA/DEC coordinates, the angle\n\
in the input catalogues must be in decimal degrees.\n", MYNAME, MYNAME);
         }
         exit(EXIT_FAILURE);
      }
      /* dump a default configuration file ---------------------------- */
      if(!strcmp(argv[i],"-d")){
         printf("# Configuration file for %s\n",MYNAME);
         printf("# ---------------------------------------------------------- #\n");
         printf("# Input catalogues                                           #\n");
         printf("# ---------------------------------------------------------- #\n");
         printf("data1          %s  # input data catalogue #1\n",    para->fileInName1);
         printf("cols1          %s,%s\t  # Column ids for data1\n",  para->data1Id[0],para->data1Id[1]);
         printf("data2          %s  # input data catalogue #2\n",    para->fileInName2);
         printf("cols2          %s,%s\t  # Column ids for data2\n",  para->data2Id[0],para->data2Id[1]);
         printf("ran1           %s\t  # input random catalogue #1\n",para->fileRanName1);
         printf("rancols1       %s,%s\t  # Column ids for ran1\n",   para->ran1Id[0],para->ran1Id[1]);
         printf("ran2           %s\t  # input random catalogue #2\n",para->fileRanName2);
         printf("rancols2       %s,%s\t  # Column ids for ran2\n",   para->ran2Id[0],para->ran2Id[1]);
         printf("fits           auto\t  # input files in fits format? [auto,yes,no(ascii)].\n");
         //   printf("coord          RADEC\t  #Coordinates: [RADEC,CART,CART3D]\n");
         //printf("                    \t  #in degrees if RADEC\n");
         printf("# ---------------------------------------------------------- #\n");
         printf("# Correlation options                                        #\n");
         printf("# ---------------------------------------------------------- #\n");
         printf("corr           auto\t # Type of correlation:\n");
         printf("                   \t # [auto,cross,gglens,auto_wp,cross_wp,auto_3D,cross_3D,number]\n");
         printf("est            ls\t # Estimator [ls,nat,ham,peebles]\n");
         printf("range          %g,%g\t # Correlation range. Dimension same as \"proj\":\n", para->min, para->max);
         printf("nbins          %d\t # Number of bins\n", para->nbins);
         printf("nbins_pi       %d\t # Number of bins for pi (for wp)\n", para->nbins_pi);
         printf("log            yes\t # Logarithmic bins [yes,no]\n");
         printf("err            jackknife # Resampling method [bootstrap,jackknife,subsample]\n");
         printf("                         # or [bootstrap2D,jackknife2D]\n");
         printf("nsub           %d\t # Number of subvolumes for resampling (power of 2)\n", para->nsamples);
         printf("nsamples       %d\t # Number of samples for resampling (= nsub for jackknife)\n", para->nsub);
         printf("seed           time\t # random seed for bootstrap\n");
         printf("OA             %g\t # Open angle for approximation (value or \"no\") \n", para->OA);
         printf("calib          no\t # Calibration factor [yes,no] (for lensing). Replace e1 by 1+m or c.\n");
         printf("RR_in          no\t # file for pre-computed RR pairs (only for clustering) \n");
         printf("RR_out         no\t # file for saving RR pairs (only for clustering) \n");
         printf("# ---------------------------------------------------------- #\n");
         printf("# Cosmology (for gal-gal correlations, w(R) and xi(rp,PI))   #\n");
         printf("# ---------------------------------------------------------- #\n");
         printf("H0             %g\t # Hubble parameter\n", para->a[0]);
         printf("Omega_M        %g\t # Relative matter density\n", para->a[1]);
         printf("Omega_L        %g\t # Relative energy density (Lambda)\n", para->a[2]);
         printf("deltaz         %g\t # For gg lensing: Zsource > Zlens + deltaz\n", para->deltaz);
         printf("pi_max         %g\t # For wp: limit for pi integration\n", para->pi_max);
         printf("# ---------------------------------------------------------- #\n");
         printf("# Output options                                             #\n");
         printf("# ---------------------------------------------------------- #\n");
         printf("proj           theta\t # Axis projection [theta, como, phys]\n");
         printf("                    \t # for auto, cross: default: theta (deg)\n");
         printf("                    \t # for auto_[3d,wp], cross_[3d,wp], gglens: default: como (Mpc)\n");
         printf("out            %s\t # Output file in ascii format\n", para->fileOutName);
         printf("cov            no\t # Covariance matrix of errors [yes,no] (in \"out\".cov)\n");
         printf("xi             no\t # xi(rp, pi) for auto_wp [yes,no] (in \"out\".xi)\n");
         printf("printTree      no\t # Output tree in \"out\".data[1,2]_tree.[ascii,fits]\n");
         printf("printTreeAndExit no\t # Same as above but exits after\n");
         printf("printSamples   no\t # Ouput samples in \"out\".samples\n");
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
            if(getStrings(line, item," ",&Ncol))
            setPara(getCharValue(item,1),getCharValue(item,2),para);
         }
         fclose(configFile);
      }
   }
   if(noconfigFile){
      if(para->verbose) fprintf(stderr,"\n%s: **ERROR** config file not found. Exiting...\n", MYNAME);
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }

   /* ----------------------------------------------------------------------
    * STEP 3: third loop over arguments. Read the command line options
    * (overwrite config file option). */
   for(i=0;i<argc;i++){
      if(*argv[i] == '-') setPara(argv[i]+1, argv[i+1], para);
   }

   /* ----------------------------------------------------------------------
    * STEP 4: readjust parameters if needed */

   /* Random seed */
   if(para->rank == MASTER && para->seed > 0){
      srand((unsigned int)para->seed);
      /*
      Debugging
      printf("setting new seed...%ld.", para->seed);
      printf(" ran = %d\n", rand());
      */
   }

   /* Default projection for gg lensing and wp */
   if( (para->corr == GGLENS || para->corr == AUTO_WP || para->corr == CROSS_WP) && para->proj == THETA ){
      para->proj = COMO;
   }

	if (para->fits == -1){
		if ( checkFileExt(para->fileInName1, ".fits") || checkFileExt(para->fileInName1, ".fit")){
			para->fits = 1;
		}else{
			para->fits = 0;
		}
	}

   /* Adjust number of dimensions */
   if(para->proj == COMO || para->proj == PHYS || para->corr == AUTO_3D || para->corr == CROSS_3D) NDIM  = 3;

   /* For number, GGLENS, or wp(rp), set NDIM = 3 and only bootsrap in 2D */
   if(para->corr == NUMBER || para->corr == AUTO_WP || para->corr == CROSS_WP){
      NDIM  = 3;
      para->resample2D = 1;
   }

   /* for 3D computation, no projection is used (-> proj=THETA) */
   if((para->proj == COMO || para->proj == PHYS) && (para->corr == AUTO_3D || para->corr == CROSS_3D)){
      para->proj = THETA;
      if(para->verbose) fprintf(stderr,"\n%s: **WARNING** proj set to theta for 3D correlation.\n", MYNAME);
   }

   /* If jackknife nsub != nsamples and check that nsamples >= nsub otherwise */
   if(para->err == JACKKNIFE && para->nsamples != para->nsub ){
      para->nsamples = para->nsub;
      if(para->verbose) fprintf(stderr,"\n%s: **WARNING** nsamples set to nsub for jackknife.\n", MYNAME);
   }
   if(para->nsamples < para->nsub ){
      if(para->verbose) fprintf(stderr,"\n%s: **WARNING** nsamples is < nsub.\n", MYNAME);
   }

   /* For GGLENS only bootsrap in 2D */
   if(para->corr == GGLENS ){
      para->resample2D = 1;
   }

   /* For 3D cartesian coordinates, angular separation replaced by 3D distance */
   if(para->corr == AUTO_3D || para->corr == CROSS_3D){
      para->distAng   = &dist3D;
   }

   /* wp(rp) integrated along pi in linear scale */
   // DEBUGGING para->Delta_pi = (para->max - para->min)/(double)para->nbins;
   para->Delta_pi = (para->pi_max - 0.0)/(double)para->nbins_pi;

   // TODO: check that all columns are defined (including for randoms)

   if(para->log){
      if(para->min < EPS){
         if(para->verbose) fprintf(stderr,"\n%s: **ERROR** logarithmic binning: min range should be > 0. Exiting...\n", MYNAME);
         MPI_Finalize();
         exit(EXIT_FAILURE);
      }
      para->min = log(para->min);
      para->max = log(para->max);
   }
   para->Delta = (para->max - para->min)/(double)para->nbins;

   if(para->err == SUBSAMPLE && para->corr != NUMBER){
      if(para->verbose) fprintf(stderr,"\n%s: **ERROR** \"subsample\" error option only defined for \"number\". Exiting...\n", MYNAME);
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }


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
		for(j=0;j<Ncol;j++) strcpy(para->data1Id[j], getCharValue(list,j+1));
		// for(j=0;j<Ncol;j++)  para->data1Id[j] = getIntValue(list,j+1);
   }else if(!strcmp(field,"data2")){
      checkArg(field,arg,para);
      strcpy(para->fileInName2,arg);
   }else if(!strcmp(field,"cols2")){
      checkArg(field,arg,para);
      getStrings(arg,list,",",&Ncol);
		for(j=0;j<Ncol;j++) strcpy(para->data2Id[j], getCharValue(list,j+1));
      // for(j=0;j<Ncol;j++)  para->data2Id[j] = getIntValue(list,j+1);
   }else if(!strcmp(field,"ran1")){
      checkArg(field,arg,para);
      strcpy(para->fileRanName1,arg);
   }else if(!strcmp(field,"rancols1")){
      checkArg(field,arg,para);
      getStrings(arg,list,",",&Ncol);
		for(j=0;j<Ncol;j++) strcpy(para->ran1Id[j], getCharValue(list,j+1));
      // for(j=0;j<Ncol;j++)  para->ran1Id[j] = getIntValue(list,j+1);
   }else if(!strcmp(field,"ran2")){
      checkArg(field,arg,para);
      strcpy(para->fileRanName2,arg);
   }else if(!strcmp(field,"rancols2")){
      checkArg(field,arg,para);
      getStrings(arg,list,",",&Ncol);
		for(j=0;j<Ncol;j++) strcpy(para->ran2Id[j], getCharValue(list,j+1));
      // for(j=0;j<Ncol;j++)  para->ran2Id[j] = getIntValue(list,j+1);
      /*
   }else if(!strcmp(field,"coord")) {
   checkArg(field,arg,para);
   if(!strcmp(arg,"RADEC"))        para->coordType = RADEC;
   else if(!strcmp(arg,"CART"))    para->coordType = CART;
   else if(!strcmp(arg,"CART3D"))  para->coordType = CART3D;
   else checkArg(field, NULL, para);
   */
}else if(!strcmp(field,"fits")) {
   checkArg(field,arg,para);
   if(!strcmp(arg,"auto"))     para->fits = -1;
   else if(!strcmp(arg,"yes")) para->fits = 1;
   else if(!strcmp(arg,"no"))  para->fits = 0;
   else checkArg(field, NULL, para);
}else if(!strcmp(field,"corr")) {
   checkArg(field, arg, para);
   if(!strcmp(arg,"cross"))         para->corr = CROSS;
   else if(!strcmp(arg,"auto"))     para->corr = AUTO;
   else if(!strcmp(arg,"gglens"))   para->corr = GGLENS;
   else if(!strcmp(arg,"auto_3D"))  para->corr = AUTO_3D;
   else if(!strcmp(arg,"cross_3D")) para->corr = CROSS_3D;
   else if(!strcmp(arg,"auto_wp"))  para->corr = AUTO_WP;
   else if(!strcmp(arg,"cross_wp")) para->corr = CROSS_WP;
   else if(!strcmp(arg,"number"))   para->corr = NUMBER;
   else checkArg(field, NULL, para);
}else if(!strcmp(field,"proj")) {
   checkArg(field,arg,para);
   if(!strcmp(arg,"theta"))      para->proj = THETA;
   else if(!strcmp(arg,"como"))  para->proj = COMO;
   else if(!strcmp(arg,"phys"))  para->proj = PHYS;
   else checkArg(field, NULL, para);
}else if(!strcmp(field,"est")) {
   if(!strcmp(arg,"ls"))           para->estimator = LS;
   else if(!strcmp(arg,"nat"))     para->estimator = NAT;
   else if(!strcmp(arg,"ham"))     para->estimator = HAM;
   else if(!strcmp(arg,"peebles")) para->estimator = PEEBLES;
   else checkArg(field, NULL, para);
}else if(!strcmp(field,"range")){
   // checkArg(field,arg,para);
   getStrings(arg,list,",",&Ncol);
   para->min = getDoubleValue(list,1);
   para->max = getDoubleValue(list,2);
}else if(!strcmp(field,"nbins")){
   checkArg(field,arg,para);
   para->nbins = atoi(arg);
}else if(!strcmp(field,"nbins_pi")){
   checkArg(field,arg,para);
   para->nbins_pi = atoi(arg);
}else if(!strcmp(field,"log")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->log = 1;
   else para->log = 0;
}else if(!strcmp(field,"err")) {
   checkArg(field,arg,para);
   if(!strcmp(arg,"jackknife"))      para->err = JACKKNIFE;
   else if(!strcmp(arg,"bootstrap")) para->err = BOOTSTRAP;
   else if(!strcmp(arg,"subsample")) para->err = SUBSAMPLE;
   else if(!strcmp(arg,"bootstrap2D")){
      para->err = BOOTSTRAP;
      para->resample2D = 1;
   } else if(!strcmp(arg,"jackknife2D")){
      para->err = JACKKNIFE;
      para->resample2D = 1;
   } else checkArg(field, NULL, para);
}else if(!strcmp(field,"nsamples")){
   checkArg(field,arg,para);
   para->nsamples = atoi(arg);
}else if(!strcmp(field,"seed")){
   checkArg(field,arg,para);
   if(strcmp(arg,"time")){
      para->seed = atoi(arg);
   }
}else if(!strcmp(field,"nsub")){
   checkArg(field,arg,para);
   para->nsub = atoi(arg);
}else if(!strcmp(field,"OA")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"no")) para->OA = -1.0;
   else para->OA     = atof(arg);
}else if(!strcmp(field,"H0")){
   checkArg(field,arg,para);
   para->a[0]   = atof(arg);
}else if(!strcmp(field,"calib")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->calib = 1;
   else para->calib = 0;
}else if(!strcmp(field,"Omega_M")){
   checkArg(field,arg,para);
   para->a[1]   = atof(arg);
}else if(!strcmp(field,"Omega_L")){
   checkArg(field,arg,para);
   para->a[2]   = atof(arg);
}else if(!strcmp(field,"deltaz")){
   checkArg(field,arg,para);
   para->deltaz  = atof(arg);
}else if(!strcmp(field,"pi_max")){
   checkArg(field,arg,para);
   para->pi_max  = atof(arg);
}else if(!strcmp(field,"weighted")){
   para->weighted = 1;
}else if(!strcmp(field,"out") || !strcmp(field,"o")){
   checkArg(field,arg,para);
   strcpy(para->fileOutName,arg);
}else if(!strcmp(field,"cov")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->cov_mat = 1;
   else para->cov_mat = 0;
}else if(!strcmp(field,"xi")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->xi = 1;
   else para->xi = 0;
}else if(!strcmp(field,"printTree")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->printTree = 1;
   else para->printTree = 0;
}else if(!strcmp(field,"printTreeAndExit")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->printTreeAndExit = 1;
   else para->printTreeAndExit = 0;
}else if(!strcmp(field,"printSamples")){
   checkArg(field,arg,para);
   if(!strcmp(arg,"yes")) para->printSamples = 1;
   else para->printSamples = 0;
}else if(!strcmp(field,"RR_in")){
   checkArg(field,arg,para);
   if(strcmp(arg,"no")){
      strcpy(para->RRInFileName, arg);
   }
}else if(!strcmp(field,"RR_out")){
   checkArg(field,arg,para);
   if(strcmp(arg,"no")){
      strcpy(para->RROutFileName, arg);
   }
}else if(!strcmp(field,"c")){
   /* do nothing, this is the option for the config file */
}else{
   /*
   if(para->verbose) fprintf(stderr,"%s: **ERROR** %s is not a valid option. Exiting...\n", MYNAME, field);
   MPI_Finalize();
   exit(EXIT_FAILURE);
   */
}

return;
}

void checkArg(char *field, char *arg, Config *para){
   /* Checks if the option has a valid argument. */
   if(arg == NULL || *arg == '-'){
      //if(arg == NULL){

      if(para->verbose) fprintf(stderr,"%s: **ERROR** option %s has no valid argument. Exiting...\n", MYNAME, field);
      MPI_Finalize();
      exit(EXIT_FAILURE);
   }
   return;
}
