#include <cmath>
#include <iostream>
#include <fstream> // file
#include <cstring> // memcpy
#include <cfloat>  // DBL_MAX
#include <vector>
#include <ctime>   // clock
#include <algorithm>
#include <vector>
#include <string>
#include <random>

int compute_num_2dpts(char *vmask,int ncams,int n3Dpts);
double compute_median(double *vec,int len);
void compute_uv(double *camera,double *point,double *uv);
void compute_right_R(double *camera,double *point,double *right_R_u,double *right_R_v);
void compute_dfdR(double*p,double *imgpts,char *vmask,double*weight,double mu,int ncams,int n3Dpts,int k,double *dfdR);
void compute_right_X(double *camera,double *point,double *right_X_u,double *right_X_v);
void compute_dfdX(double*p,double *imgpts,char *vmask,double *weight,double mu,int ncams,int n3Dpts,int k,double *dfdX);
double compute_norm1(double *a,int len);
void compute_error_vector_L2(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts,double *error_vec);
void compute_gradient(double *p,double *imgpts, double mu,char *vmask,double *weight, int ncams, int n3Dpts, double *grad);
double f(double *p,double mu,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts);
double f_L2(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts);
double Armijo_f(double *p,double mu,double* grad,double sigma,double rho,double *imgpts,char * vmask,double *weight,int ncams,int n3Dpts);
void update_p(double *p,double mu,double *grad,double sigma,double rho,double *imgpts,char *vmask,double *weight,int ncams,int n3Dpts);
void compute_error_vector(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts,double *error_vec);
void sba_L1_smoothing_grad(double *p, double *imgpts, char *vmask,double *weight, int ncams, int n3Dpts,double *newp,int step_madn);
void sba_L1_smoothing_grad(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp);
double compute_MADN(double *error_vec,int len);
void compute_madn_weight(double *error_vec,double madn_k_weight,int ncams,int n3Dpts,double Toone_std ,double *weight);
double compute_error_without0(double *p,double *weight,double *imgpts,char *vmask,int ncams,int n3Dpts);
void compute_real_2dpoints(double *p,char *vmask,int ncams,int n3Dpts,double *real_2dpoints);
void Toone(double *one_avg,double *one_std,double *mot,double *imgpts,char *vmask,int ncams,int n3Dpts);
void add_noise_gauss(double *data,double ratio_noise,int len,double *data_with_noise);
void Run(double *motstr, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp);
void compute_error_gt_all(double *newp,double *motstr,char *vmask,double *real_2dpts,double mean_one,double std_one,int ncams,int n3Dpts,double *cams_only,double *cams_3dpts,double *cams_only_recover,double *cams_3dpts_recover,int index);
int coumpute_index_outer(int len_data,double ratio_outer,int *index_outer);
void add_outer(double *data,double ratio_outer_noise,int len,int num_outer,int *index_outer,double *data_with_outer);
void random_set_cams_3dpts(int ncams,int n3Dpts,double *min_xyz,double* max_xyz,double *mot_random);
void add_noise_gauss_2dpts(double *data,double ratio_noise,int n3Dpts,double *data_with_noise);
void Toone_recover(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *motstr_recover,double *imgpts_recover);
void Toone_new(double *one_avg,double *one_std,double *mot,double *imgpts,char *vmask,int ncams,int n3Dpts);
void Toone_recover_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover);
void Toone_recover_new_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover);
void sba_L1_smoothing_grad_Wolfe(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp);
void sba_L1_smoothing_grad_Goldstein(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp);

void double_copy(double *changed_arr,double *notchanged_arr,int len);
void write_txt(int *a,int len,std::string filename);
void write_txt(double *a,int len,std::string filename);
void write_txt(char *a,int len,std::string filename);
//template<typename T> void write_txt_app(T *a,int len,std::string filename);

