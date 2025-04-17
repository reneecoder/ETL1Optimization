#include "sba_L1_smoothing.h"

void double_copy(double *changed_arr,double *notchanged_arr,int len){
    for(int i=0;i<len;i++) changed_arr[i] = notchanged_arr[i];
}

void write_txt(int *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename);
    for(i=0;i<len;i++) out<<a[i]<<std::endl;
    out.close();
}
double min_arr(double *arr,int len){
    int i;
    double out = arr[0];
    for(i=1;i<len;i++){
        if(out<arr[i]) out =arr[i];
    }
    return out;
}

double max_arr(double *arr,int len){
    int i;
    double out = arr[0];
    for(i=1;i<len;i++){
        if(out>arr[i]) out =arr[i];
    }
    return out;
}

double compute_avg(double *data,int len){
    int i;
    double out=0;
    for(i=0;i<len;i++) out+=data[i];
    out/=len;
    return out;
}

int compute_num_2dpts(char *vmask,int ncams,int n3Dpts){
    int i,num=0;
    for(i=0;i<ncams*n3Dpts;i++){
        if((int)vmask[i] == 1) num++;
    }
    return num*2;
}

void write_txt(double *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename);
    for(i=0;i<len;i++) out<<a[i]<<std::endl;
    out.close();
}

void write_txt(char *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename);
    for(i=0;i<len;i++) out<<a[i]+0<<std::endl;
    out.close();
}

void compute_mean_var(double *data,int len,double *mean,double *var){
    int i;
    double temp = 0,temp1 = 0;
    for(i=0;i<len;i++){
        temp+=data[i];
    }
    temp=temp/len;
    *mean=temp;

    for(i=0;i<len;i++){
        temp1+=(data[i]-temp)*(data[i]-temp);
    }
    temp1=temp1/len;
    *var=temp1;
}

void add_noise_gauss(double *data,double ratio_noise,int len,double *data_with_noise){
    int i;                                                                          
    double *noise,mean,var;
    std::random_device rd;

    compute_mean_var(data,len,&mean,&var);
    std::cout<<"noise mean is "<<mean<<std::endl;
    std::cout<<"noise var is "<<var<<std::endl;

    noise = new double[len];

    //定义随机数序列生成器rng，传入无符号整型作为种子，rd()会返回一个随机值（或许比time(0)好）
    std::default_random_engine rng {rd()};
    std::normal_distribution<double> dist {0.000, var*ratio_noise*0.01};

    //返回10的概率最高
    for(i=0;i<len;i++){
        noise[i] = dist(rng);
        //std::cout<<noise[i]<<std::endl;
        data_with_noise[i]=data[i] + noise[i];
    }

    delete(noise);
}

void compute_min_max(double* data,int len,double *min,double*max){
    int i;
    *min = data[0];
    *max = data[0];
    for(i=1;i<len;i++){
        if(data[i]>*max) *max = data[i];
        if(data[i]<*min) *min = data[i];
    }
}

void add_noise_gauss_2dpts(double *data,double ratio_noise,int n2Dpts,double *data_with_noise){
    int i,j;                                                                          
    double *noise,mean,var,*temp_data,min,max;
    std::random_device rd;

    temp_data = new double[n2Dpts/2];

    for(i=0;i<2;i++){
        for(j=0;j<n2Dpts/2;j++){
            temp_data[j] = data[j*2+i];
        }
        //compute_mean_var(temp_data,n3Dpts,&mean,&var);
        compute_min_max(temp_data,n2Dpts/2,&min,&max);
        std::default_random_engine rng {rd()};
        std::normal_distribution<double> dist {0.000, (1024)*ratio_noise*0.01};
        std::cout<<"max - min is "<<max-min<<std::endl;
        for(j=0;j<n2Dpts/2;j++){
            data_with_noise[j*2+i] = data[j*2+i]+dist(rng);
        }
    }

}

double compute_median(double *vec,int len){
    double output;
    std::vector<double> vector(len);
    for(int i=0;i<len;i++) vector[i]=vec[i];
    if(len % 2 ==1){
        std::nth_element(vector.begin(),vector.begin()+(len-1)/2,vector.end());
        output = vector[(len-1)/2];
    }else{
        double temp;
        std::nth_element(vector.begin(),vector.begin()+(len)/2,vector.end());
        temp = vector[(len)/2];
        std::nth_element(vector.begin(),vector.begin()+(len)/2-1,vector.end());
        output = (temp + vector[(len)/2-1])/2;
    }
    return output;
}



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

void compute_right_R(double *camera,double *point,double *right_R_u,double *right_R_v){
    //被compute_dfdR调用，因为公式太长了
    double s,alpha,beta,gamma,t0,t1,x,y,z,sin_alpha,cos_alpha,sin_beta,cos_beta,sin_gamma,cos_gamma,temp1,temp2;

    s = camera[0];
    alpha = camera[1];
    beta = camera[2];
    gamma = camera[3];
    t0 = camera[4];
    t1 = camera[5];

    x=point[0];
    y=point[1];
    z=point[2];

    sin_alpha=sin(alpha);
    cos_alpha=cos(alpha);
    sin_beta=sin(beta);
    cos_beta=cos(beta);
    sin_gamma=sin(gamma);
    cos_gamma=cos(gamma);

    temp1=(cos_beta*x+sin_alpha*sin_beta*y-cos_alpha*sin_beta*z)/(s*s);
    temp2=(cos_alpha*y+sin_alpha*z)/(s*s);
    right_R_u[0]=-cos_gamma*temp1+sin_gamma*temp2;
    right_R_v[0]=-sin_gamma*temp1-cos_gamma*temp2;

    temp1=(cos_alpha*sin_beta*y+sin_alpha*sin_beta*z)/s;
    temp2=(-sin_alpha*y+cos_alpha*z)/s;
    right_R_u[1]=cos_gamma*temp1-sin_gamma*temp2;
    right_R_v[1]=sin_gamma*temp1+cos_gamma*temp2;

    temp1=(-sin_beta*x+sin_alpha*cos_beta*y-cos_alpha*cos_beta*z)/s;
    right_R_u[2]=cos_gamma*temp1;
    right_R_v[2]=sin_gamma*temp1;

    temp1=(cos_beta*x+sin_alpha*sin_beta*y-cos_alpha*sin_beta*z)/s-t0;
    temp2=(cos_alpha*y+sin_alpha*z)/s-t1;
    right_R_u[3]= - sin_gamma*temp1-cos_gamma*temp2;
    right_R_v[3]= cos_gamma*temp1-sin_gamma*temp2;

    right_R_u[4]= - cos_gamma;
    right_R_u[5]= sin_gamma;

    right_R_v[4]= - sin_gamma;
    right_R_v[5]= - cos_gamma;
}

void compute_dfdR(double*p,double *imgpts,char *vmask,double*weight,double mu,int ncams,int n3Dpts,int k,double *dfdR){
    //计算近似函数对第k个相机参数的偏导，dfdR为输出，长度为6
    int i,j,index_temp,index_u,index_v,index_p,num_mask0;
    double temp_u,temp_v,*camera,*point3D,*uv,*right_R_u,*right_R_v;

    camera=new double[6];//当前相机参数
    point3D=new double[3];
    uv=new double[2];
    right_R_u=new double[6];
    right_R_v=new double[6];

    for(i=0;i<6;i++){
        dfdR[i]=0.0;
        camera[i] = p[k*6+i];
    }
    for(i=0;i<n3Dpts;i++){
        if((int)vmask[i*ncams+k] == 1){
            //取出第i个点的坐标
            point3D[0]=p[ncams*6+3*i];
            point3D[1]=p[ncams*6+3*i+1];
            point3D[2]=p[ncams*6+3*i+2];

            //计算当前投影
            compute_uv(camera,point3D,uv);

            //取出第i个点在第k个相机投影的坐标
            num_mask0 = 0;
            for(j=0;j<i*ncams+k;j++){
                num_mask0+=(int)vmask[j];
            }
            num_mask0 = i*ncams+k - num_mask0;
            index_u=2*i*ncams+2*k-2*num_mask0;
            index_v=2*i*ncams+2*k+1-2*num_mask0;

            //计算左项
            temp_u = weight[i*ncams+k]*(uv[0]-imgpts[index_u])/sqrt(pow((uv[0]-imgpts[index_u]),2)+mu*mu);
            temp_v = weight[i*ncams+k]*(uv[1]-imgpts[index_v])/sqrt(pow((uv[1]-imgpts[index_v]),2)+mu*mu);

            //计算右项
            compute_right_R(camera,point3D,right_R_u,right_R_v);

            //相乘后相加
            dfdR[0]+=temp_u*right_R_u[0]+temp_v*right_R_v[0];
            dfdR[1]+=temp_u*right_R_u[1]+temp_v*right_R_v[1];
            dfdR[2]+=temp_u*right_R_u[2]+temp_v*right_R_v[2];
            dfdR[3]+=temp_u*right_R_u[3]+temp_v*right_R_v[3];
            dfdR[4]+=temp_u*right_R_u[4]+temp_v*right_R_v[4];
            dfdR[5]+=temp_u*right_R_u[5]+temp_v*right_R_v[5];

            //std::cout<<temp_u<<" "<<right_R_u[0]<<" ";
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);
    delete(right_R_u);
    delete(right_R_v);
}

void compute_right_X(double *camera,double *point,double *right_X_u,double *right_X_v){
    //被compute_dfdX调用，因为公式太长了
    double s,alpha,beta,gamma,t0,t1,x,y,z,sin_alpha,cos_alpha,sin_beta,cos_beta,sin_gamma,cos_gamma,temp1,temp2;

    s = camera[0];
    alpha = camera[1];
    beta = camera[2];
    gamma = camera[3];
    t0 = camera[4];
    t1 = camera[5];

    x=point[0];
    y=point[1];
    z=point[2];

    sin_alpha=sin(alpha);
    cos_alpha=cos(alpha);
    sin_beta=sin(beta);
    cos_beta=cos(beta);
    sin_gamma=sin(gamma);
    cos_gamma=cos(gamma);

    right_X_u[0]=cos_gamma*cos_beta/s;
    right_X_u[1]=cos_gamma*sin_alpha*sin_beta/s-sin_gamma*cos_alpha/s;
    right_X_u[2]=-cos_gamma*cos_alpha*sin_beta/s-sin_gamma*sin_alpha/s;

    right_X_v[0]=sin_gamma*cos_beta/s;
    right_X_v[1]=sin_gamma*sin_alpha*sin_beta/s+cos_gamma*cos_alpha/s;
    right_X_v[2]=-sin_gamma*cos_alpha*sin_beta/s+cos_gamma*sin_alpha/s;
}

void compute_dfdX(double*p,double *imgpts,char *vmask,double *weight,double mu,int ncams,int n3Dpts,int k,double *dfdX){
    //计算近似函数对第k个三维点的偏导，dfdX为输出，长度为3
    int i,j,index_temp,index_u,index_v,index_p,num_mask0;
    double temp_u,temp_v,*camera,*point3D,*uv,*right_X_u,*right_X_v;

    camera=new double[6];//当前相机参数
    point3D=new double[3];
    uv=new double[2];
    right_X_u=new double[3];
    right_X_v=new double[3];

    for(i=0;i<3;i++){
        dfdX[i]=0;
        point3D[i]=p[6*ncams+3*k+i];
    }

    for(i=0;i<ncams;i++){

        if((int)vmask[k*ncams+i] == 1){

            //取出第i个相机的参数
            for(j=0;j<6;j++)camera[j]=p[i*6+j];

            //计算当前投影
            compute_uv(camera,point3D,uv);

            //取出第k个点在第i个相机投影的坐标
            num_mask0 = 0;
            for(j=0;j<k*ncams+i;j++){
                num_mask0+=(int)vmask[j];
            }
            num_mask0 = k*ncams+i - num_mask0;
            index_u=2*k*ncams+2*i-2*num_mask0;
            index_v=2*k*ncams+2*i+1-2*num_mask0;


            //计算左项
            temp_u = weight[k*ncams+i]* (uv[0]-imgpts[index_u])/sqrt(pow((uv[0]-imgpts[index_u]),2)+mu*mu);
            temp_v = weight[k*ncams+i]* (uv[1]-imgpts[index_v])/sqrt(pow((uv[1]-imgpts[index_v]),2)+mu*mu);

            //计算右项
            compute_right_X(camera,point3D,right_X_u,right_X_v);

            dfdX[0]+=temp_u*right_X_u[0]+temp_v*right_X_v[0];
            dfdX[1]+=temp_u*right_X_u[1]+temp_v*right_X_v[1];
            dfdX[2]+=temp_u*right_X_u[2]+temp_v*right_X_v[2];
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);
    delete(right_X_u);
    delete(right_X_v);
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

void compute_gradient(double *p,double *imgpts, double mu,char *vmask,double *weight, int ncams, int n3Dpts, double *grad) {
    /*********
    计算当前光滑近似函数的梯度
    Input：
    p为自变量，先是相机参数（s,omiga,beta,gamma,t0,t1），后是维点坐标（x，y，z）
    imgpts为点链，排列为第i个点在第j个相机下的坐标
    mu为光滑近似参数，vmask为掩码
    ncams，n3Dpts为相机个数与三维点个数
    Output：
    grad：当前近似函数的梯度
    ***********/
    int i,j,ncol;
    double *cameras,*points,*dfdR,*dfdX,du,dv;

    ncol = ncams*6+n3Dpts*3;
    memset(grad,0,ncol*sizeof(double));

    /*cameras = new double[ncams*6];
    points = new double[n3Dpts*3];
    memcpy(cameras, p,sizeof(double) * (ncams *6));
    memcpy(points, p + (ncams *6), sizeof(double) * (n3Dpts*3));*/

    dfdR = new double[6];
    memset(dfdR,0,6*sizeof(double));
    dfdX = new double[3];
    memset(dfdR,0,3*sizeof(double));

    for(i=0;i<ncams;i++){

        compute_dfdR(p,imgpts,vmask,weight,mu,ncams,n3Dpts,i,dfdR);
        memcpy(grad+i*6,dfdR,6*sizeof(double));
        /*****************
        //for(j=0;j<6;j++) grad[i*6+j]=dfdR[j];
        //for(j=0;j<6;j++) std::cout<<dfdR[j]<<" ";
        ******************/
    }

    for(i=0;i<n3Dpts;i++){
        compute_dfdX(p,imgpts,vmask,weight,mu,ncams,n3Dpts,i,dfdX);
        memcpy(grad+ncams*6+i*3,dfdX,3*sizeof(double));

    }
 
    //delete(cameras);
    //delete(points);
    delete(dfdR);
    delete(dfdX);
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

double f_L2(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts){
    //计算L2残差
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
                u_fid = imgpts[index_uv];
                index_uv++;
                v_fid = imgpts[index_uv];
                index_uv++;

                //f_out +=abs(uv[0]-u_fid)+abs(uv[1]-v_fid);
                f_out +=weight[i*ncams+j]*pow(uv[0]-u_fid,2);
                f_out +=weight[i*ncams+j]*pow(uv[1]-v_fid,2);
            }
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);

    return sqrt(f_out);
}

double Armijo_f(double *p,double mu,double* grad,double sigma,double rho,double *imgpts,char * vmask,double *weight,int ncams,int n3Dpts){
    /*************
    回退法寻找符合Armijo准则的步长
    ************/
    int i,ncols;
    bool is_continue;
    double alpha,*p_new,*p_new_temp,dot,right,f_right;

    ncols=6*ncams+3*n3Dpts;
    p_new = new double[ncols];

    alpha = 1.0;
    dot = 0.0;
    for(i=0;i<ncols;i++){
        p_new[i]=p[i]-alpha*grad[i];
        dot += (-grad[i]*grad[i]);
    }

    std::clock_t start, end;
    start = clock();//debug;

    f_right=f(p,mu,imgpts,vmask,weight,ncams,n3Dpts);
    right =f_right  + sigma*alpha*dot;
    is_continue = (f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts) > right);
    
    while( is_continue ){
        
        end=clock();
        //std::cout<<(double(end-start)/CLOCKS_PER_SEC>60)<<std::endl;
        //if(double(end-start)/CLOCKS_PER_SEC>10) break;//debug


        alpha *= rho;
        for(i=0;i<ncols;i++){
            p_new[i]=p[i]-alpha*grad[i];
        }
        right =f_right  + sigma*alpha*dot;
        is_continue = (f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts) > right);
    }

    delete(p_new);

    return alpha;
}

void update_p(double *p,double mu,double *grad,double sigma,double rho,double *imgpts,char *vmask,double *weight,int ncams,int n3Dpts){
    /*******************
    梯度法迭代，方向为负梯度，步长由Armijo准则得出
    ********************/
    int i,ncols;
    double alpha;

    ncols=6*ncams+3*n3Dpts;
    alpha = Armijo_f(p,mu,grad,sigma,rho,imgpts,vmask,weight,ncams,n3Dpts);

    for(i=0;i<ncols;i++){
        p[i]=p[i]-alpha*grad[i];
    }

}

void compute_error_vector(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts,double *error_vec){
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
                if(temp<0) temp=-temp;
                error_vec[(i*ncams+j)]=weight[i*ncams+j]*temp;
                temp=uv[1]-v_fid;
                if(temp<0) temp=-temp;
                error_vec[(i*ncams+j)]+=weight[i*ncams+j]*temp;
            }
        }
    }
    //for(i=0;i<ncams*n3Dpts;i++) std::cout<<error_vec[i]<<std::endl;
    delete(camera);
    delete(point3D);
    delete(uv);

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

void sba_L1_smoothing_grad(double *p, double *imgpts, char *vmask,double *weight, int ncams, int n3Dpts,double *newp,int step_madn) {
    int i,k,maxstep,ncols;
    double sigma,rho,mu,epsilon,omiga,norm1_grad,*grad,error_L1,*error_vec,error_L2,f_now;
    std::vector<double> error_everystep_L1,error_everystep_L2,norm1_grad_everystep;
    std::ofstream fffout;
    std::clock_t start, end;


    sigma=0.05;
    rho=0.5;
    mu=1;
    epsilon=0.001;
    omiga=1000;
    maxstep=1000;
    //手动设置的参数

    k=0;
    ncols = ncams*6 + n3Dpts*3;
    grad = new double[ncols];
    error_vec = new double[ncams*n3Dpts];

    memcpy(newp,p,sizeof(double)*ncols);

    //初始化

    //for(int i=0;i<ncols;i++) std::cout<<" "<<newp[i];

    compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
    //for(int i=0;i<ncols;i++) std::cout<<grad[i]<<" ";
    norm1_grad=compute_norm1(grad,ncols);
    //std::cout<<norm1_grad<<std::endl;
    
    //compute_error_vector(p,imgpts,vmask,weight,ncams,n3Dpts,error_vec);
    // fffout.open("out/error_everypoint_intial.txt");
    // for(i=0;i<ncams*n3Dpts;i++){
    //     fffout<<error_vec[i]<<std::endl;
    // }
    // fffout.close();
    error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);

    std::cout<<"/****************/MADN_"<<step_madn<<"_intial"<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

    start = clock();//debug;
    
    if(step_madn == 0){
        fffout.open("out/output.txt");
    }else{
        fffout.open("out/output.txt",std::ios_base::app);
    }
    
    while(norm1_grad>epsilon){
        
        //*************out every step*****************//
        end=clock();
        if(double(end-start)/CLOCKS_PER_SEC>600) break;

        error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);
        error_everystep_L1.push_back(error_L1);
        error_everystep_L2.push_back(error_L2);
        norm1_grad_everystep.push_back(norm1_grad);

        fffout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        
        //*************out every step*****************//

        update_p(newp,mu,grad,sigma,rho,imgpts,vmask,weight,ncams,n3Dpts);
        compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
        norm1_grad=compute_norm1(grad,ncols);
        if(norm1_grad<omiga*mu) mu*=sigma;
        k++;
        
        if(k>maxstep) break;
        
    }
    fffout.close();
    // compute_error_vector(newp,imgpts,vmask,weight,ncams,n3Dpts,error_vec);
    // fffout.open("out/error_everypoint_output.txt");
    // for(i=0;i<ncams*n3Dpts;i++){
    //     fffout<<error_vec[i]<<std::endl;
    // }
    // fffout.close();
    
    // fffout.open("out/error_everystep_L1.txt");
    // for(i=0;i<error_everystep_L1.size();i++){
    //     fffout<<error_everystep_L1[i]<<std::endl;
    // }
    // fffout.close();
    
    // fffout.open("out/error_everystep_L2.txt");
    // for(i=0;i<error_everystep_L2.size();i++){
    //     fffout<<error_everystep_L2[i]<<std::endl;
    // }
    // fffout.close();
    // fffout.open("out/norm_grad_everystep_L1.txt");
    // for(i=0;i<norm1_grad_everystep.size();i++){
    //     fffout<<norm1_grad_everystep[i]<<std::endl;
    // }
    // fffout.close();

    std::cout<<"/****************/MADN_"<<step_madn<<"_final"<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts)<<",error_L1 is "<<f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",error_L2 is "<<f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",norm_grad is "<<norm1_grad<<std::endl;

    delete(grad);
    delete(error_vec);
}

void sba_L1_smoothing_grad(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp) {
    int i,k,maxstep,ncols;
    double sigma,rho,mu,epsilon,omiga,norm1_grad,*grad,error_L1,*error_vec,error_L2,f_now,*weight;
    std::vector<double> error_everystep_L1,error_everystep_L2,norm1_grad_everystep;
    std::ofstream fffout,ffout1;
    std::clock_t start, end;


    sigma=0.05;
    rho=0.5;
    mu=1;
    epsilon=0.001;
    omiga=1000;
    maxstep=1000;
    //手动设置的参数

    k=0;
    ncols = ncams*6 + n3Dpts*3;
    grad = new double[ncols];
    error_vec = new double[ncams*n3Dpts];
    weight = new double[ncams*n3Dpts];
    memcpy(newp,p,sizeof(double)*ncols);
    for(i=0;i<ncams*n3Dpts;i++) weight[i] = 1.0;

    //初始化
    compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
    norm1_grad=compute_norm1(grad,ncols);

    error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);

    std::cout<<"/****************"<<"intial"<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

    start = clock();//debug;
    fffout.open("out/output.txt");
    ffout1.open("out/output_x_everystep.txt");
    
    while(norm1_grad>epsilon){
        
        //*************out every step*****************//
        end=clock();
        if(double(end-start)/CLOCKS_PER_SEC>600) break;

        
        error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);
        error_everystep_L1.push_back(error_L1);
        error_everystep_L2.push_back(error_L2);
        norm1_grad_everystep.push_back(norm1_grad);
        
        fffout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        ffout1<<"/****************/"<<k<<"/****************/"<<std::endl;
        for(i=0;i<ncams*6+n3Dpts*3;i++)ffout1<<newp[i]<<std::endl;
        //std::cout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        
        //*************out every step*****************//

        update_p(newp,mu,grad,sigma,rho,imgpts,vmask,weight,ncams,n3Dpts);
        compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
        norm1_grad=compute_norm1(grad,ncols);
        if(norm1_grad<omiga*mu) mu*=sigma;
        k++;
        
        if(k>maxstep) break;
        
    }
 
    fffout.close();
    ffout1.close();
    std::cout<<"/****************final****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts)<<",error_L1 is "<<f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",error_L2 is "<<f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",norm_grad is "<<norm1_grad<<std::endl;

    delete(grad);
    delete(error_vec);
    delete(weight);
}

double compute_MADN(double *error_vec,int len){
    int i;
    double MADN,*error_vec_temp;
    error_vec_temp=new double[len];
    MADN=compute_median(error_vec,len);
    std::cout<<"MADN temp is "<<MADN<<std::endl;
    for (i=0;i<len;i++){
        error_vec_temp[i]=(error_vec[i]-MADN)/0.6745;
        if(error_vec_temp[i]<0) error_vec_temp[i]= - error_vec_temp[i];
    }
    //for (i=0;i<len;i++) std::cout<<error_vec_temp[i]<<" ";
    MADN=compute_median(error_vec_temp,len);

    delete(error_vec_temp);
    return MADN;
}

void compute_madn_weight(double *error_vec,double madn_k_weight,int ncams,int n3Dpts,double Toone_std ,double *weight){
    /**********根据公式计算madn权重***************/
    int i,j;
    double madn_k,MADN,*error_same_cam,*error_same_3Dpts,*m_a,*m_p,m_ai,m_pj,dij;

    //for(i=0;i<ncams*n3Dpts;i++)std::cout<<error_vec[i]<<"  ";

    error_same_cam=new double[n3Dpts];
    error_same_3Dpts= new double[ncams];
    m_p=new double[ncams];
    m_a=new double[n3Dpts];

    //compute m_a and m_p
    for(i=0;i<ncams;i++){
        for(j=0;j<n3Dpts;j++){
            error_same_cam[j]=error_vec[j*ncams+i]*Toone_std;
        }
        m_p[i]=compute_median(error_same_cam,n3Dpts);
    }
    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            error_same_3Dpts[j]=error_vec[i*ncams+j]/m_p[j]*Toone_std;
        }
        m_a[i]=compute_median(error_same_3Dpts,ncams);
        if(std::isnan(m_a[i]) ) m_a[i]=0;
    }

    madn_k = 4.685*madn_k_weight;
    MADN=compute_MADN(error_vec,ncams*n3Dpts);
    std::cout<<"MADNM is "<<MADN<<std::endl;
    std::ofstream  ddebug;
    ddebug.open("dij.txt");
    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            dij=(error_vec[i*ncams+j]*Toone_std/m_p[j]-m_a[i])/(madn_k*MADN);
            if(std::isnan(dij)) dij=1.0;
            //std::cout<<dij<<"  ";
            ddebug<<error_vec[i]<<std::endl;
            weight[i*ncams+j]=(1-dij*dij)*(1-dij*dij);
            if(dij<0) weight[i*ncams+j]=1;
            if(dij>1) weight[i*ncams+j]=0;
            ddebug<<"err is "<<error_vec[i*ncams+j]<<", mpj is "<<m_p[j]<<", mai is "<<m_a[i]<<", dij is "<<dij<<std::endl;

        }
    }
    ddebug.close();

    delete(error_same_cam);
    delete(error_same_3Dpts);
    delete(m_p);
    delete(m_a);
}

double compute_error_without0(double *p,double *weight,double *imgpts,char *vmask,int ncams,int n3Dpts){
    int num_zero,i;
    double out,*weight_temp;

    weight_temp=new double[ncams*n3Dpts];
    num_zero=0;
    for(i=0;i<ncams*n3Dpts;i++){
        weight_temp[i]=1.0;
        if(weight[i]==0) {
            weight_temp[i]=0.0;
            num_zero++;
        }
    }
    out = f(p,0,imgpts,vmask,weight_temp,ncams,n3Dpts)/(ncams*n3Dpts-num_zero)/2;
    delete(weight_temp);
    return out;
}

void compute_real_2dpoints(double *p,char *vmask,int ncams,int n3Dpts,double *real_2dpoints){
    int i,j,k,index_uv;
    double *camera,*point3D,*uv,f_out,u_fid,v_fid;

    camera = new double[6];
    point3D = new double[3];
    uv = new double[2];

    index_uv = 0; 
    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                memcpy(camera,p+j*6,6*sizeof(double));
                memcpy(point3D,p+6*ncams+i*3,3*sizeof(double));
                compute_uv(camera,point3D,uv);
                real_2dpoints[index_uv++]=uv[0];
                real_2dpoints[index_uv++]=uv[1];
            }
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);

}

void Toone(double *one_avg,double *one_std,double *mot,double *imgpts,char *vmask,int ncams,int n3Dpts){
    //gui yi hua 
    int i,j,num_img,num_keys;
    double mean_temp,std_temp;
    std::vector<float> points2d;

    j=0;
    for(i=0;i<n3Dpts*ncams;i++){
        if(!((int)vmask[i] == 0)){
            points2d.push_back(imgpts[j]);
            j++;
            points2d.push_back(imgpts[j]);
            j++;
        }
    }
    for(i=0;i<ncams;i++){
        points2d.push_back(mot[i*6+4]);
        points2d.push_back(mot[i*6+5]);
        //points2d.push_back(camparams[i].s);
    }
    mean_temp = 0;
    for(i=0;i<points2d.size();i++){
        mean_temp+=points2d[i];
    }
    mean_temp=mean_temp/points2d.size();
    std_temp = 0;
    for(i=0;i<points2d.size();i++){
        std_temp+=pow((points2d[i]-mean_temp),2);
    }
    std_temp=std_temp/points2d.size();
    std_temp=sqrt(std_temp);//compute mean and std

    j=0;
    for(i=0;i<n3Dpts*ncams;i++){
        if(!((int)vmask[i] == 0)){
            imgpts[j]=(imgpts[j]-mean_temp)/std_temp;
            j++;
            imgpts[j]=(imgpts[j]-mean_temp)/std_temp;
            j++;
        }
    }
    for(i=0;i<ncams;i++){
        mot[6*i+4]=(mot[6*i+4] + mean_temp)/std_temp;
        mot[6*i+5]=(mot[6*i+5] + mean_temp)/std_temp;
    }
    for(i=6*ncams;i<6*ncams+3*n3Dpts;i++){
        mot[i]=mot[i]/std_temp;
    }
    //write_txt(mot,ncams*6+n3Dpts*3,"mot_afterone.txt");
    //write_txt(imgpts,j,"point2d_afterone.txt");

    //std::cout<<"mean is "<<mean_temp<<" ,std is "<<std_temp<<std::endl;
    *one_avg=mean_temp;
    *one_std=std_temp;

    /*for(i=0;i<num_img;i++){
        num_keys = data.m_image_data[i].m_keys.size();
        for(j=0;j<num_keys;j++){
            data.m_image_data[i].m_keys[j].m_x=(data.m_image_data[i].m_keys[j].m_x-mean_temp)/std_temp;
            data.m_image_data[i].m_keys[j].m_y=(data.m_image_data[i].m_keys[j].m_y-mean_temp)/std_temp;
        }
    }
    for(i=0;i<data.NumImages();i++){
        camparams[i].t0=(camparams[i].t0-mean_temp)/std_temp;
        camparams[i].t1=(camparams[i].t1-mean_temp)/std_temp;
        
    }
    *one_avg=mean_temp;
    *one_std=std_temp;
    */
}

void Toone_new(double *one_avg,double *one_std,double *mot,double *imgpts,char *vmask,int ncams,int n3Dpts){
    //gui yi hua 
    int i,j,num_img,num_keys;
    double mean_temp,std_temp,sin_gamma_i,cos_gamma_i;
    std::vector<float> points2d;

    j=0;
    for(i=0;i<n3Dpts*ncams;i++){
        if(!((int)vmask[i] == 0)){
            points2d.push_back(imgpts[j]);
            j++;
            points2d.push_back(imgpts[j]);
            j++;
        }
    }
    for(i=0;i<ncams;i++){
        points2d.push_back(mot[i*6+4]);
        points2d.push_back(mot[i*6+5]);
        //points2d.push_back(camparams[i].s);
    }
    mean_temp = 0;
    for(i=0;i<points2d.size();i++){
        mean_temp+=points2d[i];
    }
    mean_temp=mean_temp/points2d.size();
    std_temp = 0;
    for(i=0;i<points2d.size();i++){
        std_temp+=pow((points2d[i]-mean_temp),2);
    }
    std_temp=std_temp/points2d.size();
    std_temp=sqrt(std_temp);//compute mean and std

    j=0;
    for(i=0;i<n3Dpts*ncams;i++){
        if(!((int)vmask[i] == 0)){
            imgpts[j]=(imgpts[j]-mean_temp)/std_temp;
            j++;
            imgpts[j]=(imgpts[j]-mean_temp)/std_temp;
            j++;
        }
    }
    for(i=0;i<ncams;i++){
        sin_gamma_i = sin(mot[6*i+3]);
        cos_gamma_i = cos(mot[6*i+3]);
        mot[6*i+4]=(mot[6*i+4] + mean_temp*(cos_gamma_i+sin_gamma_i))/std_temp;
        mot[6*i+5]=(mot[6*i+5] + mean_temp*(-sin_gamma_i+cos_gamma_i))/std_temp;
    }
    for(i=6*ncams;i<6*ncams+3*n3Dpts;i++){
        mot[i]=mot[i]/std_temp;
    }
    *one_avg=mean_temp;
    *one_std=std_temp;
}

void Run(double *motstr, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp){
    /****************
    Input:
    motster：自变量，相机参数和三维点坐标，先相机参数后三维点坐标
    imgpts：标定的二维点坐标
    vmask：掩码
    ncams,n3Dpts:相机个数与三维点个数
    *****************/
    int i,step_madn;
    int temp_zero,temp_one;
    double e;
    double *weight,*error_vec1,*p,madn_k_weight,*weight_one,*weight_old,*weight_temp,Toone_std,min_weight_change,weight_change;
    std::ofstream  error_vec_out,ddebug;

    madn_k_weight = 1.0;
    step_madn = 0;
    Toone_std = 1;
    min_weight_change = 0.001;
    weight = new double[ncams*n3Dpts];
    weight_one = new double[ncams*n3Dpts];
    weight_old = new double[ncams*n3Dpts];
    weight_temp = new double[ncams*n3Dpts];
    error_vec1 = new double[ncams*n3Dpts];
    p = new double[ncams*6+n3Dpts*3];
    for(i=0;i<ncams*n3Dpts;i++) weight[i]=1.0;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i]=1.0;
    for(i=0;i<ncams*6+n3Dpts*3;i++) p[i]=motstr[i];

    // compute_madn_weight(error_vec1,ncams,n3Dpts,weight);
    // error_vec_out.open("out_error_vec/before_madn.txt");
    // for(i=0;i<ncams*n3Dpts;i++) error_vec_out<<error_vec[i]<<std::endl;
    // error_vec_out.close();

    error_vec_out.open("weightchange.txt");
    ddebug.open("errorchange.txt");
    for(step_madn=0;step_madn<5;step_madn++){
        temp_zero = 0;
        temp_one = 0;
        sba_L1_smoothing_grad(p, imgpts, vmask ,weight ,ncams,n3Dpts,newp,step_madn);
        std::cout<<"error_L1 without weight is "<<f(newp,0,imgpts,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts)/2<<std::endl;
        std::cout<<"error_L1 without zero and weight is "<<compute_error_without0(newp,weight,imgpts,vmask,ncams,n3Dpts)<<std::endl;
        //for(i=0;i<ncams*6+3*n3Dpts;i++) std::cout<<newp[i]<<std::endl;
        compute_error_vector(newp,imgpts,vmask,weight_one,ncams,n3Dpts,error_vec1);
        for(i=0;i<ncams*n3Dpts;i++) ddebug<<error_vec1[i]<<std::endl;
        for(i=0;i<ncams*n3Dpts;i++) weight_old[i]=weight[i];
        compute_madn_weight(error_vec1,madn_k_weight,ncams,n3Dpts,Toone_std,weight);
        for(i=0;i<ncams*n3Dpts;i++){
            //std::cout<<weight[i]<<"  ";
            if(weight[i] == 0) temp_zero++;
            if(weight[i] == 1) temp_one++;
        } 
        std::cout<<"zero is "<<temp_zero<<" , one is "<<temp_one<<std::endl;
        
        for(i=0;i<ncams*n3Dpts;i++) weight_temp[i]=weight_old[i]-weight[i];
        weight_change = compute_norm1(weight_temp,ncams*n3Dpts);
        std::cout<<"norm of change weight is "<<weight_change<<std::endl<<"---------------------------------------------------------------------------------"<<std::endl;
        for(i=0;i<ncams*n3Dpts;i++) error_vec_out<<weight[i]<<std::endl;
        for(i=0;i<ncams*6+n3Dpts*3;i++)p[i]=newp[i];
        if(weight_change < min_weight_change){
            std::cout<<"weight change < "<<min_weight_change<<", convergent."<<std::endl;
            break;
        }
        
    }
    error_vec_out.close();
    ddebug.close();

    // error_vec_out.open("out_error_vec/weight_madn.txt");
    // for(i=0;i<ncams*n3Dpts;i++) error_vec_out<<weight[i]<<std::endl;
    // error_vec_out.close();

    // for(i=0;i<ncams*n3Dpts;i++) weight[i]=1.0;
    // compute_error_vector(newp,imgpts,vmask,weight,ncams,n3Dpts,error_vec);
    // error_vec_out.open("out_error_vec/after_madn.txt");
    // for(i=0;i<ncams*n3Dpts;i++) error_vec_out<<error_vec[i]<<std::endl;
    // error_vec_out.close();

    std::cout<<"final L1 error with weight is "<<f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<std::endl;

    delete(p);
    delete(weight);
    delete(weight_one);
    delete(weight_old);
    delete(error_vec1);
    delete(weight_temp);
}

double compute_error_gt(double *data,double *groundtruth,int len){
    //compute the L1 error of output data and ground truth
    int i;
    double out,*temp = new double[len];
    for(i=0;i<len;i++){
        temp[i]=data[i]-groundtruth[i];
    }
    out=compute_norm1(temp,len);
    delete(temp);

    return out;
}

void Toone_recover(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *motstr_recover,double *imgpts_recover){
    int i;
    for(i=0;i<ncams;i++){
        motstr_recover[i*6]=motstr[i*6];
        motstr_recover[i*6+1]=motstr[i*6+1];
        motstr_recover[i*6+2]=motstr[i*6+2];
        motstr_recover[i*6+3]=motstr[i*6+3];
        motstr_recover[i*6+4]=motstr[i*6+4]*std_data-mean_data;
        motstr_recover[i*6+5]=motstr[i*6+5]*std_data-mean_data;
    }
    for(i=6*ncams;i<ncams*6+n3Dpts*3;i++){
        motstr_recover[i]=motstr[i]*std_data;
    }
    for(i=0;i<n2Dpts;i++){
        imgpts_recover[i] = imgpts[i]*std_data+mean_data;
    }

}

void Toone_recover_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover){
    int i;
    double sin_gamma_i,cos_gamma_i;
    for(i=0;i<ncams;i++){
        sin_gamma_i = sin(gamma_old[i]);
        cos_gamma_i = cos(gamma_old[i]);
        motstr_recover[i*6]=motstr[i*6];
        motstr_recover[i*6+1]=motstr[i*6+1];
        motstr_recover[i*6+2]=motstr[i*6+2];
        motstr_recover[i*6+3]=motstr[i*6+3];
        motstr_recover[i*6+4]=motstr[i*6+4]*std_data-mean_data*(cos_gamma_i+sin_gamma_i);
        motstr_recover[i*6+5]=motstr[i*6+5]*std_data-mean_data*(-sin_gamma_i+cos_gamma_i);
    }
    for(i=6*ncams;i<ncams*6+n3Dpts*3;i++){
        motstr_recover[i]=motstr[i]*std_data;
    }
    for(i=0;i<n2Dpts;i++){
        imgpts_recover[i] = imgpts[i]*std_data+mean_data;
    }

}

void Toone_recover_new_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover){
    int i;
    double sin_gamma_i,cos_gamma_i,sin_new,cos_new;
    for(i=0;i<ncams;i++){
        sin_gamma_i = sin(gamma_old[i]);
        cos_gamma_i = cos(gamma_old[i]);
        sin_new = sin(motstr[i*6+3]);
        cos_new = cos(motstr[i*6+3]);
        motstr_recover[i*6]=motstr[i*6];
        motstr_recover[i*6+1]=motstr[i*6+1];
        motstr_recover[i*6+2]=motstr[i*6+2];
        motstr_recover[i*6+3]=motstr[i*6+3];
        motstr_recover[i*6+4]=motstr[i*6+4]*std_data-mean_data*(cos_new+sin_new);
        motstr_recover[i*6+5]=motstr[i*6+5]*std_data-mean_data*(-sin_new+cos_new);
    }
    for(i=6*ncams;i<ncams*6+n3Dpts*3;i++){
        motstr_recover[i]=motstr[i]*std_data;
    }
    for(i=0;i<n2Dpts;i++){
        imgpts_recover[i] = imgpts[i]*std_data+mean_data;
    }

}


void compute_error_gt_all(double *newp,double *motstr,char *vmask,double *real_2dpts,double mean_one,double std_one,int ncams,int n3Dpts,double *cams_only,double *cams_3dpts,double *cams_only_recover,double *cams_3dpts_recover,int index){
    //compute the L1 error of data and ground truth, 4 errors are output
    int i;
    double *temp_p,*temp;

    temp_p = new double[ncams*6+n3Dpts*3];
    temp = new double[ncams*6+n3Dpts*3];

    cams_only[index] = compute_error_gt(newp,motstr,ncams*6);
    cams_3dpts[index] = compute_error_gt(newp,motstr,ncams*6+n3Dpts*3);

    for(i=0;i<ncams*6+n3Dpts*3;i++){
        temp_p[i] = newp[i];
        temp[i] = motstr[i];
    } 
    for(i=0;i<ncams;i++){
        temp[i*6+4] = temp[i*6+4]*std_one + mean_one;
        temp_p[i*6+5] = temp_p[i*6+5]*std_one + mean_one;
    }
    for(i=6*ncams;i<ncams*6+n3Dpts*3;i++){
        temp[i] = temp[i]*std_one;
        temp_p[i] = temp_p[i]*std_one;
    }//recover from Toone

    cams_only_recover[index] = compute_error_gt(temp_p,temp,ncams*6);
    cams_3dpts_recover[index] = compute_error_gt(temp_p,temp,ncams*6+n3Dpts*3);

    delete(temp);
    delete(temp_p);
}

int coumpute_index_outer(int len_data,double ratio_outer,int *index_outer){
    int i,num_index,temp;
    num_index = len_data*ratio_outer*0.01;
    temp = len_data / num_index;
    index_outer = new int[num_index];
    for(i=0;i<num_index;i++){
        index_outer[i] = temp*i;
        std::cout<<"outer in function is "<<index_outer[i]<<" ";
    }
    return num_index;
}
void add_outer(double *data,double range,int len,int num_outer,int *index_outer,double *data_with_outer){

    int i,index_temp,temp;
    double *noise,mean,var;
    std::random_device rd;

    compute_mean_var(data,len,&mean,&var);
    std::cout<<"var of outer is "<<var<<std::endl;
    for(i=0;i<len;i++){
        data_with_outer[i] = data[i];
    }
    if(!(num_outer == 0)){

        noise = new double[len];

        //定义随机数序列生成器rng，传入无符号整型作为种子，rd()会返回一个随机值（或许比time(0)好）
        std::default_random_engine rng {rd()};
        std::normal_distribution<double> dist {0.000, 10};//改成宽度的5%
        std::uniform_real_distribution<double> dist1{range*0.01,range*0.05};
        temp = len / num_outer;

        for(i=0;i<num_outer-1;i++){
            noise[i] = dist1(rng);
            if(dist(rng)>0) noise[i] = -noise[i];
            //std::cout<<noise[i]<<std::endl;
            index_temp = i*temp;
            index_outer[i] = index_temp;
            data_with_outer[index_temp]=data[index_temp] + noise[i];
        }

        delete(noise);
    }
}
void random_set_cams_3dpts(int ncams,int n3Dpts,double *min_xyz,double* max_xyz,double *mot_random){
    int i;
    double temp,gamma_temp;

    temp = 2.0/ncams;
    std::random_device rd;
    std::default_random_engine rng {rd()};
    std::uniform_real_distribution<double> even_x(min_xyz[0], max_xyz[0]); 
    std::uniform_real_distribution<double> even_y(min_xyz[1], max_xyz[1]); 
    std::uniform_real_distribution<double> even_z(min_xyz[2], max_xyz[2]);
    std::normal_distribution<double> dist {0.000, 1.5};
    gamma_temp = dist(rng);
    std::uniform_real_distribution<double> even_gamma(gamma_temp-0.0005, gamma_temp+0.005);
    std::normal_distribution<double> alpha_noise {0.000, temp*0.05};

    std::uniform_real_distribution<double> even_t(-20, 20);
    
    //std::cout<<"temp is "<<temp<<std::endl;
    for(i=0;i<ncams;i++){
        mot_random[i*6] = 1;
        mot_random[i*6+1] = 0;
        mot_random[i*6+2] = -1 + temp*i+alpha_noise(rng);
        mot_random[i*6+3] = even_gamma(rng);
        mot_random[i*6+4] = even_t(rng);
        mot_random[i*6+5] = even_t(rng);
    }
    for(i=0;i<n3Dpts;i++){
        mot_random[ncams*6+i*3] = even_x(rng);
        mot_random[ncams*6+i*3+1] = even_y(rng); 
        mot_random[ncams*6+i*3+2] = even_z(rng);
    }
}

double Wolfe_f(double *p,double mu,double* grad,double *imgpts,char * vmask,double *weight,int ncams,int n3Dpts){
    //改编自https://blog.csdn.net/weixin_46584887/article/details/122893222
    int i,k,ncols;
    double c1,c2,alpha,a,b,*g_0,*g_k,*p_new,f_k,l,dot,f_p,dot1,dot2;

    ncols = ncams*6+n3Dpts*3;
    c1 = 0.1;
    c2 = 0.5;
    alpha = 1;
    k = 0;
    a = 0;
    b = 9999999999;

    g_0 = new double[ncols];
    g_k = new double[ncols];
    p_new = new double[ncols];


    compute_gradient(p,imgpts,mu,vmask,weight,ncams,n3Dpts,g_0);
    f_p = f(p,mu,imgpts,vmask,weight,ncams,n3Dpts);
    for(i=0;i<ncols;i++){
        p_new[i] = p[i] +alpha*(-grad[i]);
    }
    f_k = f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts);
    compute_gradient(p_new,imgpts,mu,vmask,weight,ncams,n3Dpts,g_k);
    dot = 0;
    for(i=0;i<ncols;i++) dot += g_0[i]*(-grad[i]);
    l = f(p,mu,imgpts,vmask,weight,ncams,n3Dpts) + c1*alpha*dot;

    while (1){
        if(f_k > l){
            //std::cout<<"f_k is "<<f_k<<", l is "<<l;
            k = k+1;
            b = alpha;
            alpha = (a + b) / 2;
            for(i=0;i<ncols;i++){
                p_new[i] = p[i] +alpha*(-grad[i]);
            }
            f_k = f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts); //f_k = feval(fun, x0 + alpha * dk);

            dot = 0;
            for(i=0;i<ncols;i++) dot += (grad[i]*g_0[i]);
            l = f_p - c1*alpha*dot;
            //l = feval(fun, x0) + c1 * alpha * g0' * dk;

            compute_gradient(p_new,imgpts,mu,vmask,weight,ncams,n3Dpts,g_k);//gk = feval(grid, x0 + alpha * dk);
            continue;
        }
        dot1 = 0;
        dot2 = 0;
        for(i=0;i<ncols;i++){
            dot1 += g_k[i]*(-grad[i]);
            dot2 += c2*(-grad[i])*(g_0[i]);
        }
        //std::cout<<", dot1 is "<<dot1<<", dot2 is "<<dot2<<std::endl;
        //if (gk' * dk < c2 * g0' * dk){// Wolfe准则
        if (dot1 < dot2){
            //std::cout<<", dot1 is "<<dot1<<", dot2 is "<<dot2<<std::endl;
            k = k + 1;
            a = alpha;
            alpha = 2*alpha;
            if(alpha > (a + b) / 2) alpha = (a + b) / 2;//alpha = min([2 * alpha, (a + b) / 2]);

            for(i=0;i<ncols;i++){
                p_new[i] = p[i] +alpha*(-grad[i]);
            }
            f_k = f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts);//fk = feval(fun, x0 + alpha * dk);

            dot = 0;
            for(i=0;i<ncols;i++) dot += grad[i]*g_0[i];
            l = f_p - c1*alpha*dot;//l = feval(fun, x0) + c1 * alpha * g0' * dk;

            compute_gradient(p_new,imgpts,mu,vmask,weight,ncams,n3Dpts,g_k);
            continue;
        }
        break;
    }

	/*g0 = feval(grid, x0);	% x0处的梯度值
	fk = feval(fun, x0 + alpha * dk); 	% 函数在下一个迭代点处的目标函数值
    gk = feval(grid, x0 + alpha * dk);  % 函数在下一个迭代点处的梯度值
	l = feval(fun, x0) + c1 * alpha * g0' * dk; 	% Armjio准则
	while true
	    if fk > l
            k = k + 1;
            b = alpha;
            alpha = (a + b) / 2;
            fk = feval(fun, x0 + alpha * dk);
            l = feval(fun, x0) + c1 * alpha * g0' * dk;
            gk = feval(grid, x0 + alpha * dk);
            continue;
        end
        if gk' * dk < c2 * g0' * dk % Wolfe准则
            k = k + 1;
            a = alpha;
            alpha = min([2 * alpha, (a + b) / 2]);
            fk = feval(fun, x0 + alpha * dk);
            l = feval(fun, x0) + c1 * alpha * g0' * dk;
            gk = feval(grid, x0 + alpha * dk);
            continue;
        end
        break;
	end*/

    delete[] g_0;
    delete[] g_k;
    delete[] p_new;

    //std::cout<<"alpha is "<<alpha<<std::endl;
    return alpha;

}
void update_p_Wolfe(double *p,double mu,double *grad,double *imgpts,char *vmask,double *weight,int ncams,int n3Dpts){
    /*******************
    梯度法迭代，方向为负梯度，步长由Armijo准则得出
    ********************/
    int i,ncols;
    double alpha;

    ncols=6*ncams+3*n3Dpts;
    alpha = Wolfe_f(p,mu,grad,imgpts,vmask,weight,ncams,n3Dpts);

    for(i=0;i<ncols;i++){
        p[i]=p[i]-alpha*grad[i];
    }

}

void sba_L1_smoothing_grad_Wolfe(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp) {
    int i,k,maxstep,ncols;
    double sigma,rho,mu,epsilon,omiga,norm1_grad,*grad,error_L1,*error_vec,error_L2,f_now,*weight;
    std::vector<double> error_everystep_L1,error_everystep_L2,norm1_grad_everystep;
    std::ofstream fffout,ffout1;
    std::clock_t start, end;


    sigma=0.05;
    rho=0.5;
    mu=1;
    epsilon=0.001;
    omiga=1000;
    maxstep=1000;
    //手动设置的参数

    k=0;
    ncols = ncams*6 + n3Dpts*3;
    grad = new double[ncols];
    error_vec = new double[ncams*n3Dpts];
    weight = new double[ncams*n3Dpts];
    memcpy(newp,p,sizeof(double)*ncols);
    for(i=0;i<ncams*n3Dpts;i++) weight[i] = 1.0;

    //初始化
    compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
    norm1_grad=compute_norm1(grad,ncols);

    error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);

    std::cout<<"/****************"<<"Wolfe intial"<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

    start = clock();//debug;
    fffout.open("out/output_wolfe.txt");
    ffout1.open("out/output_x_everystep_wolfe.txt");

    while(norm1_grad>epsilon){

        //*************out every step*****************//
        end=clock();
        if(double(end-start)/CLOCKS_PER_SEC>600) break;


        error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);
        error_everystep_L1.push_back(error_L1);
        error_everystep_L2.push_back(error_L2);
        norm1_grad_everystep.push_back(norm1_grad);

        fffout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        //ffout1<<"/****************/"<<k<<"/****************/"<<std::endl;
        //for(i=0;i<ncams*6+n3Dpts*3;i++)ffout1<<newp[i]<<std::endl;
        //std::cout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

        //*************out every step*****************//

        update_p_Wolfe(newp,mu,grad,imgpts,vmask,weight,ncams,n3Dpts);
        compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
        norm1_grad=compute_norm1(grad,ncols);
        if(norm1_grad<omiga*mu) mu*=sigma;
        k++;

        if(k>maxstep) break;

    }

    fffout.close();
    ffout1.close();
    std::cout<<"/****************final****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts)<<",error_L1 is "<<f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",error_L2 is "<<f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",norm_grad is "<<norm1_grad<<std::endl;

    delete(grad);
    delete(error_vec);
    delete(weight);
}
double Goldstein_f(double *p,double mu,double* grad,double *imgpts,char * vmask,double *weight,int ncams,int n3Dpts){
    //改编自https://blog.csdn.net/weixin_46584887/article/details/122870200
    int i,k,ncols;
    double c,alpha,a,b,*g_0,*g_k,*p_new,f_0,f_k,l1,l2,dot,f_p,dot1,dot2;

    ncols = ncams*6+n3Dpts*3;
    g_0 = new double[ncols];
    g_k = new double[ncols];
    p_new = new double[ncols];


    c = 0.2;//泰勒展开式补足系数，0 < c < 1/2
    alpha = 0.5;
    k = 0;
    a = 0;
    b = 9999999999;

    compute_gradient(p,imgpts,mu,vmask,weight,ncams,n3Dpts,g_k);
    for(i=0;i<ncols;i++){
        p_new[i] = p[i] +alpha*(-grad[i]);
    }
    f_k = f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts);
    f_0 = f(p,mu,imgpts,vmask,weight,ncams,n3Dpts);

    dot = 0;
    for(i=0;i<ncols;i++) dot += g_k[i]*(-g_k[i]);
    l1 = f_0 + c*alpha*dot;
    l2 = f_0 + (1-c)*alpha*dot;

    while (1){
        if(f_k > l1){
            k = k + 1;
            b = alpha;
            alpha = (a + b) / 2;

            for(i=0;i<ncols;i++) p_new[i] = p[i] +alpha*(-grad[i]);
            f_k = f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts);

            dot = 0;
            for(i=0;i<ncols;i++) dot += g_k[i]*(-g_k[i]);
            l1 = f_0 + c*alpha*dot;
            l2 = f_0 + (1-c)*alpha*dot;
            continue;
        }
        if (f_k < l2){
            k = k + 1;
            a = alpha;
            alpha = 2 * alpha;
            if(alpha > (a + b) / 2) alpha = (a + b) / 2;

            for(i=0;i<ncols;i++) p_new[i] = p[i] +alpha*(-grad[i]);
            f_k = f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts);

            dot = 0;
            for(i=0;i<ncols;i++) dot += g_k[i]*(-g_k[i]);
            l1 = f_0 + c*alpha*dot;
            l2 = f_0 + (1-c)*alpha*dot;
            continue;
        }
        break;
    }

    /*
	c = 0.3; 	% 泰勒展开式补足系数，0 < c < 1/2
	alpha = 1; 	% 初始步长为 1
	k = 0; 		% 统计迭代次数
    a = 0; b = inf; % 二分法确定 alpha 值
	gk = feval(grid, x0);	% x0处的梯度值
	fk = feval(fun, x0 + alpha * dk); 	% 函数在下一个迭代点处的目标函数值
	l1 = feval(fun, x0) + c * alpha * gk' * dk; 	% Armjio准则
    l2 = feval(fun, x0) + (1 - c) * alpha * gk' * dk; 	% Armjio准则的补全
	while true
	    if fk > l1
            k = k + 1;
            b = alpha;
            alpha = (a + b) / 2;
            fk = feval(fun, x0 + alpha * dk);
            l1 = feval(fun, x0) + c * alpha * gk' * dk;
            l2 = feval(fun, x0) + (1 - c) * alpha * gk' * dk;
            continue;
        end
        if fk < l2
            k = k + 1;
            a = alpha;
            alpha = min([2 * alpha, (a + b) / 2]);
            fk = feval(fun, x0 + alpha * dk);
            l1 = feval(fun, x0) + c * alpha * gk' * dk;
            l2 = feval(fun, x0) + (1 - c) * alpha * gk' * dk;
            continue;
        end
        break;
	end
	xk = x0 + alpha * dk;	% 下降点
	f = feval(fun, xk);	    % 下降点处函数值
end*/
    delete[] g_0;
    delete[] g_k;
    delete[] p_new;

    //std::cout<<"alpha is "<<alpha<<std::endl;
    return alpha;
}


void update_p_Goldstein(double *p,double mu,double *grad,double *imgpts,char *vmask,double *weight,int ncams,int n3Dpts){
    /*******************
    梯度法迭代，方向为负梯度，步长由Armijo准则得出
    ********************/
    int i,ncols;
    double alpha;

    ncols=6*ncams+3*n3Dpts;
    alpha = Goldstein_f(p,mu,grad,imgpts,vmask,weight,ncams,n3Dpts);

    for(i=0;i<ncols;i++){
        p[i]=p[i]-alpha*grad[i];
    }

}

void sba_L1_smoothing_grad_Goldstein(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp){
    int i,k,maxstep,ncols;
    double sigma,rho,mu,epsilon,omiga,norm1_grad,*grad,error_L1,*error_vec,error_L2,f_now,*weight;
    std::vector<double> error_everystep_L1,error_everystep_L2,norm1_grad_everystep;
    std::ofstream fffout,ffout1;
    std::clock_t start, end;


    sigma=0.05;
    rho=0.5;
    mu=1;
    epsilon=0.001;
    omiga=1000;
    maxstep=1000;
    //手动设置的参数

    k=0;
    ncols = ncams*6 + n3Dpts*3;
    grad = new double[ncols];
    error_vec = new double[ncams*n3Dpts];
    weight = new double[ncams*n3Dpts];
    memcpy(newp,p,sizeof(double)*ncols);
    for(i=0;i<ncams*n3Dpts;i++) weight[i] = 1.0;

    //初始化
    compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
    norm1_grad=compute_norm1(grad,ncols);

    error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);

    std::cout<<"/****************"<<"Goldstein intial"<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

    start = clock();//debug;
    fffout.open("out/output_Goldstein.txt");
    ffout1.open("out/output_x_everystep_Goldstein.txt");

    while(norm1_grad>epsilon){

        //*************out every step*****************//
        end=clock();
        if(double(end-start)/CLOCKS_PER_SEC>600) break;


        error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);
        error_everystep_L1.push_back(error_L1);
        error_everystep_L2.push_back(error_L2);
        norm1_grad_everystep.push_back(norm1_grad);

        fffout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        //ffout1<<"/****************/"<<k<<"/****************/"<<std::endl;
        //for(i=0;i<ncams*6+n3Dpts*3;i++)ffout1<<newp[i]<<std::endl;
        //std::cout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

        //*************out every step*****************//

        update_p_Goldstein(newp,mu,grad,imgpts,vmask,weight,ncams,n3Dpts);
        compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
        norm1_grad=compute_norm1(grad,ncols);
        if(norm1_grad<omiga*mu) mu*=sigma;
        k++;

        if(k>maxstep) break;

    }

    fffout.close();
    ffout1.close();
    std::cout<<"/****************final****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts)<<",error_L1 is "<<f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",error_L2 is "<<f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",norm_grad is "<<norm1_grad<<std::endl;

    delete(grad);
    delete(error_vec);
    delete(weight);
}
