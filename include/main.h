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
#include "correlators.h"

/*
 *    Main routines
 */

void numberCount(Config para);
void autoCorr(Config para);
void crossCorr(Config para);
void ggCorr(Config para);


#endif
