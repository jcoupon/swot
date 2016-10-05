/*
 *    init.h
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */

#ifndef INIT_H
#define INIT_H

#include "utils.h"
//#include "tree.h"

/*
 * 	Configuration
 */




void initPara(int argc, char **argv, Config *para);
void setPara(char *field, char *arg, Config *para);
void checkArg(char *field, char *arg, Config *para);


#endif
