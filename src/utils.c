#include "utils.h"

/*
 *    utils.c
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */


/*
 *		Utils - Numeric
 */

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
      while((line[i] == *delimit || line[i] == '\t') && (line[i] != '\0' || line[i] != '#' || line[i] != '\n')) i++;
      begin = i;
      while((line[i] != *delimit && line[i] != '\t') && (line[i] != '\0' && line[i] != '#' && line[i] != '\n')) i++;
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




void comment(const Config para, char *commentString){
   if(para.verbose){
      fflush(stderr);
      fprintf(stderr, "%s", commentString);
      fflush(stderr);
   }
   return;
}



/*
 * 	Point routines
 */


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
   a.w[i]    = b.w[j];
   if(para.corr == GGLENS){
      a.zerr[i] = b.zerr[j];
      a.e1[i]   = b.e1[j];
      a.e2[i]   = b.e2[j];
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
   a.w[i]    = 0.0;
   if(para.corr == GGLENS){
      a.zerr[i] = 0.0;
      a.e1[i]   = 0.0;
      a.e2[i]   = 0.0;
   }
   /* sum up */
   for(n=0;n<point.N;n++){
      for(dim=0;dim<NDIM;dim++) a.x[NDIM*i+dim] += point.x[NDIM*n+dim];
      a.w[i]    +=  point.w[n];
      if(para.corr == GGLENS){
         a.zerr[i] +=  point.zerr[n];
         a.e1[i]   +=  point.e1[n];
         a.e2[i]   +=  point.e2[n];
      }
   }
   /* divide by N */
   for(dim=0;dim<NDIM;dim++) a.x[NDIM*i+dim] /= (double)point.N;
   a.w[i]    /= (double)point.N;
   if(para.corr == GGLENS){
      a.zerr[i] /= (double)point.N;
      a.e1[i]   /= (double)point.N;
      a.e2[i]   /= (double)point.N;
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
   a->w = b.w + shift;
   if(para.corr == GGLENS){
      a->zerr = b.zerr + shift;
      a->e1   = b.e1   + shift;
      a->e2   = b.e2   + shift;
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
   tmp           = point.w[i];
   point.w[i]    = point.w[j];
   point.w[j]    = tmp;
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
   }

   return;
}


Point createPoint(const Config para, long N){
   Point point;

   point.N   = N;
   point.dim = -1;

   point.x = (double *)malloc(N*NDIM*sizeof(double));
   point.w = (double *)malloc(N*sizeof(double));
   if(para.corr == GGLENS){
      point.zerr = (double *)malloc(N*sizeof(double));
      point.e1   = (double *)malloc(N*sizeof(double));
      point.e2   = (double *)malloc(N*sizeof(double));
   }

   return point;
}

void freePoint(const Config para, Point point){

   free(point.x);
   free(point.w);
   if(para.corr == GGLENS){
      free(point.zerr);
      free(point.e1);
      free(point.e2);
   }

   return;
}

Point readCat(const Config para, char *fileInName, int id[NIDSMAX], int weighted){
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
      if(para.verbose) fprintf(stderr,"\n%s: **WARNING** the number of points is lower than the number of cpus.\n", MYNAME);
      //  exit(EXIT_FAILURE);
   }

   Point data = createPoint(para, N);

   n = 0;
   while(fgets(line,NFIELD*NCHAR,fileIn) != NULL){

      if(getStrings(line,item," ",&Ncol)){
         /* read the value at the column definied by "id[i]" */
         for(dim=0;dim<NDIM;dim++) data.x[NDIM*n+dim] = getDoubleValue(item,id[dim]);
         if(weighted){
            data.w[n] = getDoubleValue(item,id[NDIM+0]);
         }else{
            data.w[n] = 1.0;
         }
         if(para.corr == GGLENS){
            Ncol -= 4;
            //Check for positive values here ??
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





double distAngCart(const Tree *a, const long *i, const Tree *b, const long *j){
   /* Returns the distance between nodes
   * a[i] and b[j]. Cartesian coordinates in 2D.
   */

   double d0 = (b->point.x[NDIM*(*j)+0] - a->point.x[NDIM*(*i)+0]);
   double d1 = (b->point.x[NDIM*(*j)+1] - a->point.x[NDIM*(*i)+1]);

   return sqrt(d0*d0 + d1*d1);
}

double distAngPointCart(const Config *para, const Point *a, const long *i, const Point *b, const long *j){
   /* returns the angular distance between points
   * a[i] and b[i]. Cartesian coordinates in 2D.
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


double dist3D(const Tree *a, const long *i, const Tree *b, const long *j){
   /* Returns the distance between nodes
   * a[i] and b[j]. Cartesian coordinates in 3D.
   */

   double d0 = (b->point.x[NDIM*(*j)+0] - a->point.x[NDIM*(*i)+0]);
   double d1 = (b->point.x[NDIM*(*j)+1] - a->point.x[NDIM*(*i)+1]);
   double d2 = (b->point.x[NDIM*(*j)+2] - a->point.x[NDIM*(*i)+2]);

   return sqrt(d0*d0 + d1*d1 + d2*d2);
}
