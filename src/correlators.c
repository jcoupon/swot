#include "correlators.h"

/*
 *    correlators.c
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */


/*
 *		correlation function routines
 */


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
      switch(estimator){
         case LS:  /*   Landy and Szalay */
            if(D1D2.NN[para.nbins*l+i] > 0
            && D1R2.NN[para.nbins*l+i] > 0
            && D2R1.NN[para.nbins*l+i] > 0
            && R1R2.NN[para.nbins*l+i] > 0){
               result  =  Norm1*D1D2.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i];
               result += -Norm2*D1R2.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i];
               result += -Norm3*D2R1.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i] + 1.0;
            }
            break;
         case NAT: /*   Natural */
            if(D1D2.NN[para.nbins*l+i] > 0
            && R1R2.NN[para.nbins*l+i] > 0){
               result = Norm1*D1D2.NN[para.nbins*l+i]/R1R2.NN[para.nbins*l+i] - 1.0;
            }
            break;
         case HAM: /*   Hamilton */
            if(D1D2.NN[para.nbins*l+i] > 0
            && D1R2.NN[para.nbins*l+i] > 0
            && D2R1.NN[para.nbins*l+i] > 0
            && R1R2.NN[para.nbins*l+i] > 0){
               result = Norm4*D1D2.NN[para.nbins*l+i]*R1R2.NN[para.nbins*l+i]/(D1R2.NN[para.nbins*l+i]*D2R1.NN[para.nbins*l+i]) - 1.0;
            }
            break;
         case PEEBLES: /*  Peebles */
            if(D1D2.NN[para.nbins*l+i] > 0
            && D1R2.NN[para.nbins*l+i] > 0){
               result = Norm5*D1D2.NN[para.nbins*l+i]/D1R2.NN[para.nbins*l+i] - 1.0;
            }
            break;
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
   // double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/((R1R2.N2[l]-1)*(D1D2.N2[l]-1));
   // double Norm5 = R1R2.N1[l]/(D1D2.N2[l]-1.0);
   double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]))/(D1D2.N1[l]*(D1D2.N2[l]));
   double Norm2 = (R1R2.N1[l])/D1D2.N1[l];
   double Norm3 = (R1R2.N2[l])/D1D2.N2[l];
   double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/(R1R2.N2[l]*D1D2.N2[l]);
   double Norm5 = R1R2.N1[l]/D1D2.N2[l];

   double result = 0.0, sum = 0.0;

   if(i > -1){
      // if( D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0
      //&& D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0
      //&& D2R1.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0
      //&& R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0){

      if(R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0){
         switch(estimator){
            case LS:  /*   Landy and Szalay */
	            result  =  Norm1*D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)];
	            result += -Norm2*D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)];
	            result += -Norm3*D2R1.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] + 1.0;
	            break;
            case NAT: /*    Natural */
	            result = Norm1*D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] - 1.0;
	            break;
            case HAM: /*   Hamilton */
	            result = Norm4*D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)]*R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)]/(D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)]*D2R1.NN[j + para.nbins*(i + para.nbins_pi*l)]) - 1.0;
	            break;
            case PEEBLES: /*  Peebles */
	            result = Norm5*D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)]/D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] - 1.0;
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
      for(i=0; i<para.nbins_pi; i++){

         sum = 0.0;
         //     if( D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0
         //	  && D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0
         //	  && D2R1.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0
         //	  && R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0){

         if( R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] > 0.0){
            switch(estimator){
               case LS:  /* Landy and Szalay */
	               sum  =  Norm1*D1D2.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)];
	               sum += -Norm2*D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)];
	               sum += -Norm3*D2R1.NN[j + para.nbins*(i + para.nbins_pi*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] + 1.0;
	               break;
               case NAT: /* Natural */
	               sum = Norm1*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] - 1.0;
	               break;
               case HAM: /* Hamilton */
	               sum = Norm4*D1D2.NN[j + para.nbins*(i + para.nbins*l)]*R1R2.NN[j + para.nbins*(i + para.nbins_pi*l)]/(D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)]*D2R1.NN[j + para.nbins*(i + para.nbins_pi*l)]) - 1.0;
	               break;
               case PEEBLES: /* Peebles */
	               sum = Norm5*D1D2.NN[j + para.nbins*(i + para.nbins*l)]/D1R2.NN[j + para.nbins*(i + para.nbins_pi*l)] - 1.0;
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
   // double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/((R1R2.N2[l]-1)*(D1D2.N2[l]-1));
   // double Norm5 = R1R2.N1[l]/(D1D2.N2[l]-1.0);
   double Norm1 = (R1R2.N1[l]*(R1R2.N2[l]))/(D1D2.N1[l]*(D1D2.N2[l]));
   double Norm2 = (R1R2.N1[l])/D1D2.N1[l];
   double Norm3 = (R1R2.N2[l])/D1D2.N2[l];
   double Norm4 = (D1D2.N2[l]*R1R2.N2[l])/(R1R2.N2[l]*D1D2.N2[l]);
   double Norm5 = R1R2.N1[l]/D1D2.N2[l];


   double result = 0.0;
   switch(estimator){
      case LS:  /*   Landy and Szalay */
         if(D1D2.NN_s[para.nbins*l+i] > 0
         && D1R2.NN_s[para.nbins*l+i] > 0
         && D2R1.NN_s[para.nbins*l+i] > 0
         && R1R2.NN_s[para.nbins*l+i] > 0){
            result  =  Norm1*D1D2.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i];
            result += -Norm2*D1R2.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i];
            result += -Norm3*D2R1.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i] + 1.0;
         }
         break;
      case NAT: /*   Natural */
         if(D1D2.NN_s[para.nbins*l+i] > 0
         && R1R2.NN_s[para.nbins*l+i] > 0){
            result = Norm1*D1D2.NN_s[para.nbins*l+i]/R1R2.NN_s[para.nbins*l+i] - 1.0;
         }
         break;
      case HAM: /*   Hamilton */
         if(D1D2.NN_s[para.nbins*l+i] > 0
         && D1R2.NN_s[para.nbins*l+i] > 0
         && D2R1.NN_s[para.nbins*l+i] > 0
         && R1R2.NN_s[para.nbins*l+i] > 0){
            result = Norm4*D1D2.NN_s[para.nbins*l+i]*R1R2.NN_s[para.nbins*l+i]/(D1R2.NN_s[para.nbins*l+i]*D2R1.NN_s[para.nbins*l+i]) - 1.0;
         }
         break;
      case PEEBLES: /*  Peebles */
         if(D1D2.NN_s[para.nbins*l+i] > 0
         && D1R2.NN_s[para.nbins*l+i] > 0){
            result = Norm5*D1D2.NN_s[para.nbins*l+i]/D1R2.NN_s[para.nbins*l+i] - 1.0;
         }
         break;
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
      // DEBUGGING
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

Result NpairsWp(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall){
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
      result.NN = (double *)malloc(para->nbins*para->nbins_pi*(para->nsamples+1)*sizeof(double));
      for(k = 0; k < para->nbins*para->nbins_pi*(para->nsamples+1); k++) result.NN[k] = 0.0;
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
         NpairsWp(para, tree1, tree1->left[i],  tree2, tree2->left[j],  0);
         NpairsWp(para, tree1, tree1->right[i], tree2, tree2->left[j],  0);
         NpairsWp(para, tree1, tree1->left[i],  tree2, tree2->right[j], 0);
         NpairsWp(para, tree1, tree1->right[i], tree2, tree2->right[j], 0);
      }else{
         NpairsWp(para, tree1, tree1->left[i], tree2,  j,  0);
         NpairsWp(para, tree1, tree1->right[i], tree2, j,  0);
      }
   }else if(node(tree2,j) && tree2->r[j]/d > para->OA){
      NpairsWp(para, tree1, i, tree2, tree2->left[j],   0);
      NpairsWp(para, tree1, i, tree2,  tree2->right[j], 0);
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

      if(0 <= k && k < para->nbins_pi && 0 <= m && m < para->nbins){
         NN = tree1->N[i]*tree2->N[j];
         result.NN[m + para->nbins*k] += NN;
         for(l=0;l<para->nsamples;l++){
            result.NN[m + para->nbins*(k + (para->nbins_pi*(l+1)))] += NN*tree1->w[para->nsamples*i + l]*tree2->w[para->nsamples*j + l];
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
      for(k = 0; k < para->nbins*para->nbins_pi*(para->nsamples+1); k++){
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
