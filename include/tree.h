/*
 *    tree.h
 *    swot (Super W Of Theta)  mpi version
 *    Jean Coupon (2012-2016)
 */

#ifndef TREE_H
#define TREE_H

#include "utils.h"
#include "init.h"


/*
 *		Tree routines
 */

Tree buildTree(const Config *para, Point *data, Mask *mask, int dim, int firstCall);
void resample(const Config *para, const Point *data, int dim, Mask *mask, Mask const *limits, int firstCall);
void freeMask(const Config para, Mask mask);
void freeTree(Config para, Tree tree);
long splitTree(const Config *para, const Tree *tree1, const long root, const long Ncpu, int firstCall);
void printTreeFits(const Config para, char *fileOutName, const Tree tree, long i, long NLeaf, int firstCall);
void printTree(const Config para, char *fileOutName, const Tree tree, long i, long NLeaf, int firstcall);
long countNodes(long N, long NLeaf);


#endif
