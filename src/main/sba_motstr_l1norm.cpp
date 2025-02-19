#include "sba_motstr_l1norm.h"
#include <cmath>
#include <iostream>
#include <fstream> // file
#include <cstring> // memcpy
#include <cfloat>  // DBL_MAX
#include <vector>
#include <ctime>   // clock

// #include "img_projs.h"
//#include "Eigen/Dense"  //  dense matrix
//#include "Eigen/Sparse" // sparse matrix
//#include "Eigen/IterativeLinearSolvers"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RMatrixXd;

Eigen::VectorXd LinEqnsSolverLDLT(const Eigen::MatrixXd A, const Eigen::VectorXd b) {
    return A.ldlt().solve(b);
}

// Sparse QR
Eigen::VectorXd LinEqnsSolverQR(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b) {
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        // decomposition failed
    }
    return solver.solve(b);
}

// A must be square matrix
Eigen::VectorXd LinEqnsSolverLU(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd b) {
    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        // puts("decomposition failed");
        // return Eigen::VectorXd::Ones(A.cols());
    }
    return solver.solve(b);
}

void calcImgParallelProj(const double param[6], const double M[3], double n[2])
{
	double s, alpha, beta, gamma, t0, t1;
	double X, Y, Z;
	s = param[0]; alpha = param[1]; beta = param[2]; gamma = param[3]; t0 = param[4]; t1 = param[5];
	X = M[0]; Y = M[1]; Z = M[2];
	
	double cos_alpha = cos(alpha);
	double sin_alpha = sin(alpha);
	double cos_beta  = cos(beta);
	double sin_beta  = sin(beta);
	double cos_gamma = cos(gamma);
	double sin_gamma = sin(gamma);
	
	double tm1, tm2;
	tm1 = (cos_beta *X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)/s - t0;
	tm2 = (cos_alpha*Y+sin_alpha*Z)/s                               - t1;
	n[0] = cos_gamma*tm1-sin_gamma*tm2;
	n[1] = sin_gamma*tm1+cos_gamma*tm2;
}

void calcImgParallelProjJac(double param[6], double M[3], double jacmP[2][6], double jacmS[2][3]) // Ref
{
	double s, alpha, beta, gamma, t0, t1;
	double X, Y, Z;
	s = param[0]; alpha = param[1]; beta = param[2]; gamma = param[3]; t0 = param[4]; t1 = param[5];
	X = M[0]; Y = M[1]; Z = M[2];
	
	double cos_alpha = cos(alpha);
	double sin_alpha = sin(alpha);
	double cos_beta = cos(beta);
	double sin_beta = sin(beta);
	double cos_gamma = cos(gamma);
	double sin_gamma = sin(gamma);
	double s_s = 1/s;
	jacmP[0][0] = (-cos_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)+sin_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
	jacmP[0][1] = (cos_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)-sin_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
	jacmP[0][2] = cos_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
	jacmP[0][3] = -sin_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-cos_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
	jacmP[0][4] = -cos_gamma;
	jacmP[0][5] = sin_gamma;
	jacmP[1][0] = (-sin_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)-cos_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
	jacmP[1][1] = (sin_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)+cos_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
	jacmP[1][2] = sin_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
	jacmP[1][3] = cos_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-sin_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
	jacmP[1][4] = -sin_gamma;
	jacmP[1][5] = -cos_gamma;
	
	jacmS[0][0] = cos_gamma*cos_beta*s_s;
	jacmS[0][1] = (cos_gamma*sin_alpha*sin_beta-sin_gamma*cos_alpha)*s_s;
	jacmS[0][2] = (-cos_gamma*cos_alpha*sin_beta-sin_gamma*sin_alpha)*s_s;
	jacmS[1][0] = sin_gamma*cos_beta*s_s;
	jacmS[1][1] = (sin_gamma*sin_alpha*sin_beta+cos_gamma*cos_alpha)*s_s;
	jacmS[1][2] = (-sin_gamma*cos_alpha*sin_beta+cos_gamma*sin_alpha)*s_s;
}

// ok
void JacMat(double *p, char *vmask, int ncams, int n3Dpts, /* return value */double *sjacm) {
    double *images = new double[ncams*6];
    double *points = new double[n3Dpts*3];
    memcpy(images, p,              sizeof(double) * (ncams *6));
    memcpy(points, p + (ncams *6), sizeof(double) * (n3Dpts*3));
    int ncol = ncams*6+n3Dpts*3;
    double param[6], M[3], jacmP[2][6], jacmS[2][3];
    for (int i = 0; i < n3Dpts; i++) {
        for (int j = 0; j < ncams; j++) {
            if (vmask[i*ncams+j]) { // ith 3D point is visible in jth camera
                memcpy(param, images+6*j, sizeof(double)*6);
                memcpy(M,     points+3*i, sizeof(double)*3);
                // **********************************There has a question*****************************************
                //                  calccalcImgParallelProjJac -> calcImgParallelProjJac
                calcImgParallelProjJac(param, M, jacmP, jacmS);
                // ***********************************************************************************************
                for (int r = 0; r < 2; r++) {
                    for (int c = 0; c < 6; c++) {
                        // sjacm[(i*2*ncams+j*2 + r)][(j*6 + c)] = jacmP[r][c];
                        sjacm[(i*2*ncams+j*2 + r)*ncol+(j*6 + c)] = jacmP[r][c];
                    }
                }
                for (int r = 0; r < 2; r++) {
                    for (int c = 0; c < 3; c++) {
                        sjacm[(i*2*ncams+j*2 + r)*ncol+(6*ncams+i*3 + c)] = jacmS[r][c];
                    }
                }
            }
        }
    }
}

Eigen::VectorXd f(int *p, char *vmask, int ncams, int n3Dpts) {
    double *images = new double[ncams *6];
    double *points = new double[n3Dpts*3];
    memcpy(images, p,              sizeof(double) * (ncams *6));
    memcpy(points, p + (ncams *6), sizeof(double) * (n3Dpts*3));
    
    std::vector<double> f;
    double param[6], M[3], n[2];
    for (int i = 0; i < n3Dpts; i++) {
        for (int j = 0; j < ncams; j++) {
            if (vmask[i*ncams+j]) {
                memcpy(param, images+6*j, sizeof(double)*6);
                memcpy(M,     points+3*i, sizeof(double)*3);
                calcImgParallelProj(param, M, n);
                f.push_back(n[0]);
                f.push_back(n[1]);
            } else {
                f.push_back(0.0);
                f.push_back(0.0);
            }
        }
    }
    Eigen::Map<Eigen::VectorXd> vec_f(f.data(), f.size());
    return vec_f;
}


// ok
void ResErr(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts, double *r) {
    double *images = new double[ncams *6];
    double *points = new double[n3Dpts*3];
    memcpy(images, p,              sizeof(double) * (ncams *6));
    memcpy(points, p + (ncams *6), sizeof(double) * (n3Dpts*3));
    
    int gap = 0;
    std::vector<double> f, u;
    double param[6], M[3], n[2];
    for (int i = 0; i < n3Dpts; i++) {
        for (int j = 0; j < ncams; j++) {
            if (vmask[i*ncams+j]) {
                memcpy(param, images+6*j, sizeof(double)*6);
                memcpy(M,     points+3*i, sizeof(double)*3);
                calcImgParallelProj(param, M, n);
                f.push_back(n[0]);
                f.push_back(n[1]);

                u.push_back(imgpts[i*(ncams*2)+2*j + 0 - gap*2]);
                u.push_back(imgpts[i*(ncams*2)+2*j + 1 - gap*2]);
           } else {
                f.push_back(0.0);
                f.push_back(0.0);
                u.push_back(0.0);
                u.push_back(0.0);
                gap++;
            }
        }
    }
    for (int i = 0; i < ncams*2*n3Dpts; i++) r[i] = f[i]-u[i];
}

// ok
// return (double *)c1,c2,d1,d2,d,g
void ComputeC1C2D1D2(double *jac, double *r, double *x, double t, int ncams, int n3Dpts, 
/* return value: */double *c1, double *c2, double *d1, double *d2, double *d, double *g) {
    Eigen::Map<RMatrixXd> mat_jac(jac, ncams*2*n3Dpts, ncams*6+n3Dpts*3);
    Eigen::SparseMatrix<double>  sp_jac = mat_jac.sparseView();
    
    Eigen::Map<Eigen::VectorXd> vec_r (r, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_x (x, ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 
    Eigen::VectorXd dp = vec_x.head(ncams*6+n3Dpts*3); 
    Eigen::VectorXd  s = vec_x.tail(ncams*2*n3Dpts); 
    Eigen::VectorXd jac_dp_r;

    jac_dp_r = sp_jac*dp + vec_r;

    for (int i = 0; i < ncams*2*n3Dpts; i++) { 
        c1[i] = 1.0 / (s[i] - jac_dp_r[i]);
        c2[i] = 1.0 / (s[i] + jac_dp_r[i]);
        d1[i] = c1[i]*c1[i] + c2[i]*c2[i];
        d2[i] = c2[i]*c2[i] - c1[i]*c1[i];
  
        d[i] = 2.0                       /(s[i]*s[i] + jac_dp_r[i]*jac_dp_r[i]);
        g[i] = 2.0*jac_dp_r[i]*(s[i]*t-1)/(s[i]*s[i] + jac_dp_r[i]*jac_dp_r[i]);
    }
}

// ok
// return the Newton step (double *)dx and decrement (double)lambda
void ComputeNewtonStepAndDecrement(double *jac, double *r, double *x, double t, int ncams, int n3Dpts, double *dx, double &lambda) {
    Eigen::Map<RMatrixXd> mat_jac(jac, ncams*2*n3Dpts, ncams*6+n3Dpts*3);
    Eigen::SparseMatrix<double>  sp_jac   = mat_jac.sparseView();
    Eigen::SparseMatrix<double>  sp_jac_t =  sp_jac.transpose();

    double *arr_c1 = new double[ncams*2*n3Dpts];
    double *arr_c2 = new double[ncams*2*n3Dpts];
    double *arr_d1 = new double[ncams*2*n3Dpts];
    double *arr_d2 = new double[ncams*2*n3Dpts];
    double *arr_d  = new double[ncams*2*n3Dpts];
    double *arr_g  = new double[ncams*2*n3Dpts];

    double *arr_ddp = new double[ncams*6+n3Dpts*3];
    double *arr_ds  = new double[ncams*2*n3Dpts];

    Eigen::Map<Eigen::VectorXd> vec_c1(arr_c1, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_c2(arr_c2, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_d1(arr_d1, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_d2(arr_d2, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_d (arr_d,  ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_g (arr_g,  ncams*2*n3Dpts); 

    Eigen::Map<Eigen::VectorXd> vec_ddp(arr_ddp, ncams*6+n3Dpts*3); 
    Eigen::Map<Eigen::VectorXd> vec_ds (arr_ds,  ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_dx (    dx,  ncams*6+n3Dpts*3 + ncams*2*n3Dpts);
    Eigen::VectorXd            antigrad(         ncams*6+n3Dpts*3 + ncams*2*n3Dpts);

    ComputeC1C2D1D2(jac, r, x, t, ncams, n3Dpts, arr_c1, arr_c2, arr_d1, arr_d2, arr_d, arr_g);

    Eigen::SparseMatrix<double> sp_diag_d = Eigen::MatrixXd(vec_d.asDiagonal()).sparseView();
    Eigen::SparseMatrix<double> sp_a1; 
    Eigen::VectorXd                b1, b2; 

    sp_a1 =  sp_jac_t * sp_diag_d * sp_jac;
    //sp_a1 =  sp_jac * sp_diag_d * sp_jac_t;
       b1 = -sp_jac_t * vec_g;
    //sp_a1 = sp_a1.transpose();
    vec_ddp = LinEqnsSolverLDLT(Eigen::MatrixXd(sp_a1), b1); // cond(a1) too large
    
    b2 = (-t)*Eigen::VectorXd::Ones(ncams*2*n3Dpts) + vec_c1 + vec_c2 
       + (-1)*Eigen::SparseMatrix<double>(Eigen::MatrixXd(vec_d2.asDiagonal()).sparseView())*sp_jac*vec_ddp;

    for (int i = 0; i < ncams*2*n3Dpts; i++) {
        vec_ds[i] = 1.0/vec_d1[i] * b2[i];
    }
    vec_dx.head(ncams*6+n3Dpts*3) = vec_ddp;
    vec_dx.tail(ncams*2*n3Dpts)   = vec_ds;

    antigrad.head(ncams*6+n3Dpts*3) = -sp_jac_t*(vec_c1-vec_c2);
    antigrad.tail(ncams*2*n3Dpts)   = (-t)*Eigen::VectorXd::Ones(ncams*2*n3Dpts) + vec_c1 + vec_c2;

    lambda = antigrad.dot(vec_dx);
}

// ok
double LogBarrierFunc(double *jac, double *r, double *x, double t, int ncams, int n3Dpts) {
    Eigen::Map<RMatrixXd> mat_jac(jac, ncams*2*n3Dpts, ncams*6+n3Dpts*3);
    Eigen::SparseMatrix<double>  sp_jac = mat_jac.sparseView();

    Eigen::Map<Eigen::VectorXd> vec_r (r,                    ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_x (x, ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 

    Eigen::VectorXd dp = vec_x.head(ncams*6+n3Dpts*3); 
    Eigen::VectorXd  s = vec_x.tail(ncams*2*n3Dpts); 
    
    Eigen::VectorXd jac_dp;
    Eigen::VectorXd Ax(2*(ncams*2*n3Dpts));
    Eigen::VectorXd b (2*(ncams*2*n3Dpts));

    jac_dp = sp_jac*dp;
    Ax.head(ncams*2*n3Dpts) =  jac_dp - s;
    Ax.tail(ncams*2*n3Dpts) = -jac_dp - s;

    b. head(ncams*2*n3Dpts) = -vec_r;
    b. tail(ncams*2*n3Dpts) =  vec_r;
    
    double phi = t * s.sum();
    for (int i = 0; i < 2*(ncams*2*n3Dpts); i++) {
        if (b[i]-Ax[i] <= 0) return DBL_MAX;
        phi -= log(b[i]-Ax[i]);
    }
    return phi;
}

// ok
char BLSCond(double *jac, double *r, double *x, double *dx, double alpha, double step, double t, double lambda, int ncams, int n3Dpts) {
    double *nx = new double[ncams*6+n3Dpts*3 + ncams*2*n3Dpts]; // nx = x+step*dx
    Eigen::Map<Eigen::VectorXd> vec_x (x,  ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_dx(dx, ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_nx(nx, ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 
    vec_nx = vec_x + step*vec_dx;
    return LogBarrierFunc(jac, r, nx, t, ncams, n3Dpts) > LogBarrierFunc(jac, r, x, t, ncams, n3Dpts) - alpha*step*lambda;
}

// ok
double BackLineSearch(double *jac, double *r, double *x, double *dx, double t, double lambda, int ncams, int n3Dpts) {
    double alpha = 0.01, beta = 0.5, step = 1.0;
    while (BLSCond(jac, r, x, dx, alpha, step, t, lambda, ncams, n3Dpts)) step *= beta;
    return step;
}

// ok
// return (double *)x
void NewtonsMethod(double *jac, double *r, double *x, double t, int ncams, int n3Dpts) {
    double *dx = new double[ncams*6+n3Dpts*3 + ncams*2*n3Dpts], lambda;
    Eigen::Map<Eigen::VectorXd> vec_x (x,  ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_dx(dx, ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 
    double eps = 1e-5;

    double last = 0.0;

    while (1) {
        ComputeNewtonStepAndDecrement(jac, r, x, t, ncams, n3Dpts, dx, lambda);
        // printf("lambda=%f\n", lambda);
        if (fabs(lambda/2.0 - last) < 0.1*eps) {
            // puts("Newton's method failed!");
            break;
        }

        if (lambda/2.0 < eps) break;
        double step = BackLineSearch(jac, r, x, dx, t, lambda, ncams, n3Dpts);
        vec_x += step * vec_dx;
  
        last = lambda/2.0;
    }
}

// return (double *)x
void BarrierMethod(double *jac, double *r, int ncams, int n3Dpts, double *x) {
    Eigen::Map<Eigen::VectorXd> vec_r (r, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_x (x, ncams*6+n3Dpts*3 + ncams*2*n3Dpts); 

    vec_x.head(ncams*6+n3Dpts*3) =                    Eigen::VectorXd::Zero(ncams*6+n3Dpts*3);
    vec_x.tail(ncams*2*n3Dpts)   = vec_r.cwiseAbs() + Eigen::VectorXd::Ones(ncams*2*n3Dpts);

    double m = 2.0*(ncams*2*n3Dpts), t = m/vec_x.tail(ncams*2*n3Dpts).sum(), miu = 10.0, eps = 1e-5;
    // printf("m=%f,t=%f\n", m, t);
    while (1) {
        NewtonsMethod(jac, r, x, t, ncams, n3Dpts);
        if (m/t < eps) break;
        t *= miu;
    }
}

// Input p, u(calc in imgpts), eta, lambda
void sba_motstr_l1norm(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts, double *rk, int &cnt) {
    double *jac = new double[(ncams*2*n3Dpts) * (ncams*6+n3Dpts*3)];

    double *r   = new double[ncams*2*n3Dpts];
//    double *rk  = new double[ncams*2*n3Dpts];
    
    double *dp  = new double[ncams*6+n3Dpts*3];
    double *np  = new double[ncams*6+n3Dpts*3]; // stand for "next p", aka np = p+dp
    double * x  = new double[ncams*6+n3Dpts*3 + ncams*2*n3Dpts];
    
    Eigen::Map<Eigen::VectorXd> vec_r  (r,  ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_rk (rk, ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_p  (p,  ncams*6+n3Dpts*3); 
    Eigen::Map<Eigen::VectorXd> vec_dp (dp, ncams*6+n3Dpts*3);      
    Eigen::Map<Eigen::VectorXd> vec_np (np, ncams*6+n3Dpts*3);      
    Eigen::Map<Eigen::VectorXd> vec_x  (x,  ncams*6+n3Dpts*3 + ncams*2*n3Dpts);

    Eigen::Map<RMatrixXd>       mat_jac(jac, ncams*2*n3Dpts, ncams*6+n3Dpts*3);


    //***********************************************************************************************************
    // int cnt = 0;
    cnt = 0;
    double int_err, out_error, end_err_value = 1e-3;
    std::vector<double> out_err_list;
    //***********************************************************************************************************

    ResErr(p, imgpts, vmask, ncams, n3Dpts, rk); // rk=r0 for k=0

    //***********************************************************************************************************
    printf("initial error L1: %f\n", vec_rk.lpNorm<1>()/vec_rk.size());
    printf("initial error L2: %f\n", vec_rk.lpNorm<2>()/vec_rk.size());
    
    for (int i = 0; i < vec_rk.size(); i++) {
        if (vec_rk[i] == 0.0) cnt++;
    }
    // printf("count=%d\n", cnt);
    //***********************************************************************************************************

    double eta = 1e-10, lambda = 2.0;
    char converge = 0;

    while (!converge) {
        JacMat(p, vmask, ncams, n3Dpts, jac);     // Get Jacobian matrix at p
        BarrierMethod(jac, rk, ncams, n3Dpts, x);  // Get x = [dp, s]
        vec_dp = vec_x.head(ncams*6+n3Dpts*3);
        vec_np = vec_p + vec_dp;
        ResErr(np, imgpts, vmask, ncams, n3Dpts, r);
        while (vec_r.lpNorm<1>() >= vec_rk.lpNorm<1>()) {
            vec_dp /= lambda;
            if (vec_dp.lpNorm<1>() < eta) {
                converge = 1;
                break;
            }
            vec_np = vec_p + vec_dp;
            ResErr(np, imgpts, vmask, ncams, n3Dpts, r);
            //*******************************************************************************************************
            int_err = vec_r.lpNorm<1>()/(vec_r.size()-cnt);
            if (int_err < end_err_value) {
                vec_p  = vec_p + vec_dp;
                vec_rk = vec_r;
                return;
            } else {
                printf("inter error: %f\n", int_err);
            }
            //*******************************************************************************************************
        }
        vec_p  = vec_p + vec_dp;
        vec_rk = vec_r;
        //*******************************************************************************************************
        out_error = vec_rk.lpNorm<1>() / (vec_rk.size() - cnt);
        out_err_list.push_back(out_error);
        // cnt_outerr++;
        int len_err_list = out_err_list.size();
        if (len_err_list >= 10) {
            double *out_err_list_tail = new double[10];
            for (int i = 0; i < 10; i++) out_err_list_tail[i] = out_err_list[len_err_list - 10 + i];
            double max, min;
            max = min = out_err_list_tail[0];
            for (int i = 0; i < 10; i++) {
                if (max <= out_err_list_tail[i]) max = out_err_list_tail[i];
                if (min >= out_err_list_tail[i]) min = out_err_list_tail[i];
            }
            if ( (max - min) <= 1e-1) {
                vec_p  = vec_p + vec_dp;
                vec_rk = vec_r;
                return;
            }
        } else {
            double max, min;
            max = min = out_err_list[0];
            for (int i = 0; i < len_err_list; i++) {
                if (max <= out_err_list[i]) max = out_err_list[i];
                if (min >= out_err_list[i]) min = out_err_list[i];
            }
            if ( (max - min) <= 1e-2) {
                vec_p  = vec_p + vec_dp;
                vec_rk = vec_r;
                return;
            }
        }
        if (out_error < end_err_value) {
            vec_p  = vec_p + vec_dp;
            vec_rk = vec_r;
            return;
        } else {
            printf("out error: %f\n", out_error);
        }
        //***********************************************************************************************************
    }
}
double compute_norm11(double *a,int len){
    int i;
    double output;
    output=0;
    for(i=0;i<len;i++){
        if(a[i]>0)output+=a[i];
        if(a[i]<0)output=output-a[i];
        //output+=abs(a[i]);
    }
    return output;
}


double Run_inner(double *motstr, double *imgpts, char *vmask, int ncams, int n3Dpts) {
    /* motstr  ==  parameter_init_P0  in python code
       imgpts  ==  list_m             in python code */
    // **********************************TEST***********************************************
    // **********************************DONE***********************************************
    int cnt = 0;
    double *r = new double[ncams*2*n3Dpts];
    Eigen::Map<Eigen::VectorXd> vec_r     (r,      ncams*2*n3Dpts); 
    Eigen::Map<Eigen::VectorXd> vec_motstr(motstr, ncams*6+n3Dpts*3);
    
    sba_motstr_l1norm(motstr, imgpts, vmask, ncams, n3Dpts, r, cnt);
    
    /*printf("r0_count: %d\n", cnt);
    printf("len_r_final: %d\n", vec_r.size());
    puts("r_final:");
    std::cout << vec_r.transpose() << std::endl;
    puts("final_camera_L1_new:");
    std::cout << vec_motstr.head(ncams *6).transpose() << std::endl;
    puts("final_point_L1_new:");
    std::cout << vec_motstr.tail(n3Dpts*3).transpose() << std::endl;*/
    //printf("final error L1 of inner: %f\n", vec_r.lpNorm<1>()/vec_r.size());
    // printf("final error L2: %f\n", vec_r.lpNorm<2>()/vec_r.size());

    // std::ofstream  fout;
    // fout.open("afterinner.txt");
    // for(int i=0;i<ncams*6+n3Dpts*3;i++){
    //     fout<<vec_motstr[i]<< std::endl;
    // }
    std::cout<<"norm1 is  "<<compute_norm11(r,ncams*2*n3Dpts)<<std::endl;
    return vec_r.lpNorm<1>()/vec_r.size();
    // fout.close();

}
/*
int main() {
    std::ifstream  fin;
    
    int ncams = 64, n3Dpts = 37;
    double *motstr = new double[ncams*6+n3Dpts*3];
    double *imgpts = new double[n3Dpts*ncams*2];
    char   *vmask  = new char  [n3Dpts*ncams];
    memset(imgpts, 0, sizeof(double) * (ncams*2*n3Dpts));
    
    fin.open("camparams_ba_beforeLM.txt", std::ios::in);
    for (int i = 0; i < ncams*6+n3Dpts*3; i++) fin >> motstr[i];
    fin.close();
//    for (int i = 0; i < ncams*6+n3Dpts*3; i++) std::cout << motstr[i] << std::endl;

    fin.open("projections_ba_beforeLM.txt", std::ios::in);
    for (int i = 0; i < 2176*2; i++) fin >> imgpts[i];
    fin.close();
//     for (int i = 0; i < 2112*2; i++) std::cout << imgpts[i] << std::endl;

    fin.open("vmask_ba_beforeLM.txt", std::ios::in);
    for (int i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
    }
    fin.close();
//    for (int i = 0; i < n3Dpts*ncams; i++) std::cout << vmask[i] << std::endl;

    std::clock_t start, end;
    start = clock(); 
    Run(motstr, imgpts, vmask, ncams, n3Dpts);
    end = clock();
    std::cout << "time = " << double(end-start)/CLOCKS_PER_SEC << "s" << std::endl;  
    return 0;
}*/
