/*
 *    correlators.h
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */

#ifndef CORRELATORS_H
#define CORRELATORS_H

#include "utils.h"
#include "tree.h"
#include "init.h"



/*
 *		correlation function routines
 */


double wTheta(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int l);
double wp(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int j, int l);
double xis(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int l);

Result Nobjects(const Config *para, const Tree *tree1, const long i, int firstCall);
Result Npairs(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall);
Result NpairsWp(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall);
Result gg(const Config *para,const Tree *lens, const long i, const Tree *source, const long j, int firstCall);
void corrLensSource(const Config *para, const Tree *lens, long i, const Tree *source, long j, double deltaTheta, Result result);



#endif
