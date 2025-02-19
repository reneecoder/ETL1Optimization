#ifndef SBA_MOTSTR_L1NORM_H // if not define
#define SBA_MOTSTR_L1NORM_H

#include <vector>
#include "Eigen/Dense"  //  dense matrix
#include "Eigen/Sparse" // sparse matrix
#include "Eigen/IterativeLinearSolvers"



double Run_inner(double *motstr, double *imgpts, char *vmask, int ncams, int n3Dpts);

#endif