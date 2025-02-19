//
// Created by 徐子赫 on 2023/4/4.
//

#include "BA.h"
template<>
bool ReprojectionError<6>::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{
    util::point2d p;
    CameraParameters cam;
    cam.s=parameters[0][0];
    cam.alpha=parameters[0][1];
    cam.beta=parameters[0][2];
    cam.gamma=parameters[0][3];
    cam.t0=parameters[0][4];
    cam.t1=parameters[0][5];

    double t0 = cam.t0;
    double t1 = cam.t1;

    Eigen::Map<const Eigen::Vector3d> point(parameters[1]);

    double X, Y, Z;
    X=parameters[1][0];
    Y=parameters[1][1];
    Z=parameters[1][2];

    p=cam.Project(point);

    residuals[0] = p.x  - _observation_x;
    residuals[1] = p.y  - _observation_y;

    double cos_alpha = cos(cam.alpha);
    double sin_alpha = sin(cam.alpha);
    double cos_beta = cos(cam.beta);
    double sin_beta = sin(cam.beta);
    double cos_gamma = cos(cam.gamma);
    double sin_gamma = sin(cam.gamma);
    double s_s = 1/cam.s;

    if(jacobians !=NULL)
    {
        if(jacobians[0] != NULL)
        {
//            Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor> > J_se3(jacobians[0]);
            Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> J(jacobians[0]);
            J(0,0)=(-cos_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)+sin_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;;
            J(0,1) = (cos_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)-sin_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
            J(0,2) = cos_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
            J(0,3) = -sin_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-cos_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
            J(0,4) = -cos_gamma;
            J(0,5) = sin_gamma;
            J(1,0) = (-sin_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)-cos_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
            J(1,1) = (sin_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)+cos_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
            J(1,2) = sin_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
            J(1,3) = cos_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-sin_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
            J(1,4) = -sin_gamma;
            J(1,5) = -cos_gamma;
        }

        if(jacobians[1] != NULL)
        {
            jacobians[1][0] = cos_gamma*cos_beta*s_s;
            jacobians[1][1] = (cos_gamma*sin_alpha*sin_beta-sin_gamma*cos_alpha)*s_s;
            jacobians[1][2] = (-cos_gamma*cos_alpha*sin_beta-sin_gamma*sin_alpha)*s_s;
            jacobians[1][3] = sin_gamma*cos_beta*s_s;
            jacobians[1][4] = (sin_gamma*sin_alpha*sin_beta+cos_gamma*cos_alpha)*s_s;
            jacobians[1][5] = (-sin_gamma*cos_alpha*sin_beta+cos_gamma*sin_alpha)*s_s;
        }
    }

    return true;

}

util::point2d CameraParameters::Project(const Eigen::Vector3d& Point3D) {

    double cos_alpha;
    double sin_alpha;
    double cos_beta;
    double sin_beta;
    double cos_gamma;
    double sin_gamma;

    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);
    cos_beta = cos(beta);
    sin_beta = sin(beta);
    cos_gamma = cos(gamma);
    sin_gamma = sin(gamma);

    double a1,a2;
    double u,v;

    a1 = (cos_beta * Point3D[0] + sin_alpha * sin_beta * Point3D[1] - cos_alpha * sin_beta * Point3D[2]) / s - t0;
    a2 = (cos_alpha * Point3D[1] + sin_alpha * Point3D[2]) / s - t1;

    u = cos_gamma * a1 - sin_gamma * a2;
    v = sin_gamma * a1 + cos_gamma * a2;

    util::point2d p;
    p.x=u;
    p.y=v;

    return p;
}

/*
void compute_uv(double *camera,double *point,double *uv){
    //计算投影u
    double u,v,s,alpha,beta,gamma,t0,t1,x,y,z,temp1,temp2;

    s = camera[0];
    alpha = camera[1];
    beta = camera[2];
    gamma = camera[3];
    t0 = camera[4];
    t1 = camera[5];

    x=point[0];
    y=point[1];
    z=point[2];

    temp1=(cos(beta)*x+sin(alpha)*sin(beta)*y-cos(alpha)*sin(beta)*z)/s-t0;
    temp2=(cos(alpha)*y+sin(alpha)*z)/s-t1;

    u=cos(gamma)*temp1-sin(gamma)*temp2;
    v=sin(gamma)*temp1+cos(gamma)*temp2;

    uv[0]=u;
    uv[1]=v;
}

void compute_error_vector_L2(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts,double *error_vec){
    //计算每个点的残差值
    int i,j,k,index_uv;
    double *camera,*point3D,*uv,f_out,u_fid,v_fid,temp;

    camera = new double[6];
    point3D = new double[3];
    uv = new double[2];

    f_out =0;
    index_uv = 0;
    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            error_vec[(i*ncams+j)]=0.0;
            if((int)vmask[i*ncams+j] == 1){
                memcpy(camera,p+j*6,6*sizeof(double));
                memcpy(point3D,p+6*ncams+i*3,3*sizeof(double));
                compute_uv(camera,point3D,uv);
                u_fid = imgpts[index_uv];
                index_uv++;
                v_fid = imgpts[index_uv];
                index_uv++;
                temp=uv[0]-u_fid;
                error_vec[(i*ncams+j)]=weight[i*ncams+j]*temp*temp;
                temp=uv[1]-v_fid;
                error_vec[(i*ncams+j)]+=weight[i*ncams+j]*temp*temp;
                error_vec[(i*ncams+j)] = sqrt(error_vec[(i*ncams+j)]);
            }
        }
    }
    //for(i=0;i<ncams*n3Dpts;i++) std::cout<<error_vec[i]<<std::endl;
    delete(camera);
    delete(point3D);
    delete(uv);

}

double f(double *p,double mu,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts){
    //计算光滑近似函数的值
    int i,j,k,index_uv;
    double *camera,*point3D,*uv,f_out,u_fid,v_fid;

    camera = new double[6];
    point3D = new double[3];
    uv = new double[2];

    f_out =0;
    index_uv = 0;
    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                memcpy(camera,p+j*6,6*sizeof(double));
                memcpy(point3D,p+6*ncams+i*3,3*sizeof(double));
                compute_uv(camera,point3D,uv);
                //std::cout<<" "<<uv[0]<<" "<<uv[1]<<std::endl;
                u_fid = imgpts[index_uv];
                index_uv++;
                v_fid = imgpts[index_uv];
                index_uv++;

                //f_out +=abs(uv[0]-u_fid)+abs(uv[1]-v_fid);
                f_out +=weight[i*ncams+j]*sqrt(pow(uv[0]-u_fid,2)+mu*mu);
                f_out +=weight[i*ncams+j]*sqrt(pow(uv[1]-v_fid,2)+mu*mu);
                //std::cout<<"f"<<f_out<<"u"<<uv[0]<<"v"<<uv[1]<<" ";
                //std::cout<<"f"<<f_out<<"f+"<<sqrt(pow(uv[0]-u_fid,2))<<" ";
            }
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);

    return f_out;
}

double compute_norm1(double *a,int len){
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

*/

void LM_BA(double *Camera_noise,double *imgpts_state,char* vmask,int ncams,int n3Dpts,double *newp){
    int gap,i;
    gap=0;
    ceres::Problem problem;
    double tmp=0;
    double x_o,y_o,dij;

    PosePointParametersBlock init_states;
    init_states.create(ncams, n3Dpts);

    for(i=0;i<ncams;i++)
    {
        init_states.values[i*6]=Camera_noise[i*6];
        init_states.values[i*6+1]=Camera_noise[i*6+1];
        init_states.values[i*6+2]=Camera_noise[i*6+2];
        init_states.values[i*6+3]=Camera_noise[i*6+3];
        init_states.values[i*6+4]=Camera_noise[i*6+4];
        init_states.values[i*6+5]=Camera_noise[i*6+5];
    }
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++)
    {
        init_states.values[i]=Camera_noise[i];
    }

    for(int i=0; i<n3Dpts; i++){
        for(int j=0;j<ncams; j++){
            dij=1;
            if(vmask[i*ncams+j]){
                x_o=imgpts_state[i*(ncams*2)+2*j+0-gap*2];
                y_o=imgpts_state[i*(ncams*2)+2*j+1-gap*2];
                ceres::CostFunction* cost_function;
                cost_function=new ReprojectionError<6>(x_o, y_o);
                ceres::LossFunction* lossFunc = new MADN_loss(dij);
                problem.AddResidualBlock(cost_function,lossFunc, init_states.pose(j), init_states.point(i));
            }
            else{
                gap++;
            }
        }
    }

    ceres::Solver::Options options;
    //配置增量方程的解法
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = false;
    options.max_num_iterations = 1000;

    //第三步，创建Summary对象用于输出迭代结果
    ceres::Solver::Summary summary;

//    ceres::Solver(options,&problem,&summary);
    ceres::Solve(options, &problem, &summary);
    for(i=0;i<ncams*6+n3Dpts*3;i++) newp[i] = init_states.values[i];
    //std::cout << summary.BriefReport() << "\n";

}
void LM_BA_weight(double *Camera_noise,double *imgpts_state,char* vmask,double* weight_pts2d,int ncams,int n3Dpts,double *newp){
    int gap,i;
    gap=0;
    ceres::Problem problem;
    double tmp=0;
    double x_o,y_o,dij;

    PosePointParametersBlock init_states;
    init_states.create(ncams, n3Dpts);

    for(i=0;i<ncams;i++)
    {
        init_states.values[i*6]=Camera_noise[i*6];
        init_states.values[i*6+1]=Camera_noise[i*6+1];
        init_states.values[i*6+2]=Camera_noise[i*6+2];
        init_states.values[i*6+3]=Camera_noise[i*6+3];
        init_states.values[i*6+4]=Camera_noise[i*6+4];
        init_states.values[i*6+5]=Camera_noise[i*6+5];
    }
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++)
    {
        init_states.values[i]=Camera_noise[i];
    }

    for(int i=0; i<n3Dpts; i++){
        for(int j=0;j<ncams; j++){
            dij=weight_pts2d[i*ncams+j];
            if(vmask[i*ncams+j]){
                x_o=imgpts_state[i*(ncams*2)+2*j+0-gap*2];
                y_o=imgpts_state[i*(ncams*2)+2*j+1-gap*2];
                ceres::CostFunction* cost_function;
                cost_function=new ReprojectionError<6>(x_o, y_o);
                ceres::LossFunction* lossFunc = new MADN_loss(dij);
                problem.AddResidualBlock(cost_function,lossFunc, init_states.pose(j), init_states.point(i));
            }
            else{
                gap++;
            }
        }
    }
    ceres::Solver::Options options;
    //配置增量方程的解法
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = false;
    options.max_num_iterations = 1000;

    //第三步，创建Summary对象用于输出迭代结果
    ceres::Solver::Summary summary;

//    ceres::Solver(options,&problem,&summary);
    ceres::Solve(options, &problem, &summary);
    for(i=0;i<ncams*6+n3Dpts*3;i++) newp[i] = init_states.values[i];
    //std::cout << summary.BriefReport() << "\n";

}

void BA_IRLS(double *p,double *pts2d,char* vmask,int ncams,int n3Dpts,double *newp){
    int i,j,max_step;
    double *p_temp,*error_vec,*weight_one,*weight,*p_change,epsilon;

    max_step = 500;
    epsilon = 0.001;

    p_temp = new double[ncams*6+n3Dpts*3];
    p_change = new double[ncams*6+n3Dpts*3];
    error_vec = new double[ncams*n3Dpts];
    weight = new double[ncams*n3Dpts];
    weight_one = new double[ncams*n3Dpts];
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;

    LM_BA(p,pts2d,vmask,ncams,n3Dpts,p_temp);//第一次不带权的LM
    compute_error_vector_L2(p_temp,pts2d,vmask,weight_one,ncams,n3Dpts,error_vec);
    for(i=0;i<ncams*n3Dpts;i++) {
        weight[i] = 0;
        if(!(error_vec[i] == 0))   weight[i] = 1/error_vec[i];
    }//计算第一次权重
    for(i=0;i<ncams*6+n3Dpts*3;i++) newp[i] = p_temp[i];//更新自变量

    for(i=0;i<max_step;i++){

        //std::cout<<"L1 error is "<<f(newp,0,pts2d,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts)/2<<std::endl;//输出每次迭代的L1误差

        LM_BA_weight(newp,pts2d,vmask,weight,ncams,n3Dpts,p_temp);//带权重的LM

        for(j=0;j<ncams*6+n3Dpts*3;j++) {
            p_change[j] = newp[j]-p_temp[j];
            newp[j] = p_temp[j];
        }
        if(compute_norm1(p_change,ncams*6+n3Dpts*3)<epsilon) break;//终止条件之一。自变量变化值（L1范数）小于epsilon

        compute_error_vector_L2(p_temp,pts2d,vmask,weight_one,ncams,n3Dpts,error_vec);
        for(j=0;j<ncams*n3Dpts;j++) {
            weight[j] = 0;
            if(!(error_vec[j] == 0))   weight[j] = 1/error_vec[j];
        }//计算下一次的权重

    }

    delete[] p_temp;
    delete[] error_vec;
    delete[] weight_one;
    delete[] weight;
}
