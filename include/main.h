/*
 *    main.h
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */

#ifndef MAIN_H
#define MAIN_H

#include "utils.h"
#include "tree.h"
#include "init.h"

/*
 *    Main routines
 */

void numberCount(Config para);
void autoCorr(Config para);
void crossCorr(Config para);
void ggCorr(Config para);

double wTheta(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int l);
double wp(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int j, int l);
double xis(const Config para, int estimator, Result D1D2, Result R1R2, Result D1R1, Result D2R2, int i, int l);

Result Nobjects(const Config *para, const Tree *tree1, const long i, int firstCall);
Result Npairs(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall);
Result Npairs3D(const Config *para, const Tree *tree1, const long i, const Tree *tree2, const long j, int firstCall);
Result gg(const Config *para,const Tree *lens, const long i, const Tree *source, const long j, int firstCall);
void corrLensSource(const Config *para, const Tree *lens, long i, const Tree *source, long j, double deltaTheta, Result result);

#endif
