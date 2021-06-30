#include <stdio.h>
#include "../include/matrix.h"
#include "../include/kalman.h"
#include <math.h>
#include <random>
#include <fstream>

using namespace std;



/*C matrix evolution function on iter_param iteration*/
TMatrix C_k(TMatrix state_params, nat iter_param)

{

    nat n = 3,m=1;
    double sample_rate = 0.01;
  
    TMatrix Ck = create_matrix(m,n);
    double A = get_data(state_params,0,0);
    double xi = get_data(state_params,1,0);
    double w = get_data(state_params,2,0);
    double t = iter_param*sample_rate;
 
    double exp_term = exp(-xi*w*t);
    double sin_term = sin(w*t+atan(1/xi));
    double cos_term = cos(w*t+atan(1/xi));

    insert_data(Ck,0,0,1-exp_term*sin_term);
    insert_data(Ck,0,1,t*A*w*exp_term*sin_term + t*A*exp_term*cos_term/(1+pow(xi,2)));
    insert_data(Ck,0,2,t*A*xi*exp_term*sin_term - t*A*exp_term*cos_term);
  
    return Ck;
  
}

/*state evolution function*/
TMatrix x_k(TMatrix state_params, nat iter_param)

{
    nat n = 3;
    TMatrix xk = create_matrix(n,1);
  
    /*restriction*/
    double xi = get_data(state_params,1,0);
    if (xi < 0){
      insert_data(state_params,1,0,0.0001);
    }
    if (xi > 1){
      insert_data(state_params,1,0,0.9999);
    }
  
    /*updating state*/
    for (nat i = 0 ; i < n ; i++){
      insert_data(xk,i,0,get_data(state_params,i,0));
    }
    return xk;
  
}


/*output observation evolution function on iter param iteration*/
TMatrix y_k(TMatrix state_params, nat iter_param)

{

    nat m = 1;
    double sample_rate = 0.01;

    TMatrix yk = create_matrix(m,m);
    double A = get_data(state_params,0,0);
    double xi = get_data(state_params,1,0);
    double w = get_data(state_params,2,0);
    double t = iter_param*sample_rate;

    double exp_term = exp(-xi*w*t);
    double sin_term = sin(w*t+atan(1/xi));

    insert_data(yk,0,0,A-A*exp_term*sin_term);

    return yk;
  
}


/*Phi evolution function on iter_param iteration (identity matrix)*/
TMatrix Phi_k(TMatrix state_params, nat iter_param)

{
    nat n = 3;
    TMatrix Phi = identity(n);
    return Phi;
}


/*G evolution function on iter_param iteration (zero matrix)*/
TMatrix G_k(TMatrix state_params, nat iter_param)

{
    nat n = 3;
    TMatrix G = create_matrix(n,n);
    return G;
}


/*U evolution function on iter_param iteration (identity matrix)*/
TMatrix U_k(TMatrix state_params, nat iter_param)

{
    nat m = 1;
    TMatrix U = identity(m);
    return U;
}


int main()


{

    /*n = states size, m = output size*/
    /*n_samples = number of iterations of EKF*/
    nat n = 3, m = 1;
    nat n_samples = 200;

    /*States prior*/
    TMatrix x_prior = create_matrix(n,1);
    double xk_data[n*1] = {0.01,0.01,0.01};
    fill_matrix(x_prior,xk_data);

    /*P prior*/
    TMatrix P_prior = create_matrix(n,n);
    double Pp_data[9] = {100,0,0,0,1,0,0,0,100};
    fill_matrix(P_prior,Pp_data);

    /*Q, defualt zero matrix*/
    TMatrix Q = create_matrix(n,n);

    /*R*/
    TMatrix R = create_matrix(m,m);
    insert_data(R,0,0,4);

    /*I*/
    TMatrix I = identity(n);
   
    /*Matrix that will contain true signal for comparation*/
    TMatrix signal = create_matrix(m,n_samples);
    
    /*y*/
    TMatrix y = create_matrix(m,n_samples);
    double y_data[n_samples];
    double y_signal[n_samples];

    double A = 8;
    double xi = 0.1;
    double w = 15;
    double sample_rate = 0.01;
    double exp_term,sin_term;
    double xi_sqrt = sqrt(1-pow(xi,2));

    /*noise*/
    const double mean = 0.0;
    const double stddev = 2;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean,stddev);

    /*Signal generation*/
    for (nat i = 0; i < n_samples;i++){


      exp_term = exp(-xi*w*i*sample_rate)/xi_sqrt;
      sin_term = sin(w * xi_sqrt *  i*sample_rate + atan(xi_sqrt/xi));


      y_signal[i] = A - A*exp_term*sin_term;
      y_data[i] = A-A*exp_term*sin_term+dist(generator);

    }

    fill_matrix(y,y_data);
    fill_matrix(signal,y_signal);


    /*creating EKF_func*/
    TEKF_func yk = create_kalman_f(x_prior,0,y_k);
    TEKF_func xk = create_kalman_f(x_prior,0,x_k);
    TEKF_func G = create_kalman_f(x_prior,0,G_k);
    TEKF_func C = create_kalman_f(x_prior,0,C_k);
    TEKF_func Phi = create_kalman_f(x_prior,0,Phi_k);
    TEKF_func U = create_kalman_f(x_prior,0,U_k);


    /*EKF*/
    TMatrix out_matrix = EKF(I,P_prior,Q,R,y,x_prior,C,Phi,G,U,xk,yk);

    
    /*estimated signal and states extraction from EKF output*/
    TMatrix estimated_signal = get_signal(out_matrix,1);
    TMatrix states = get_state(out_matrix,1);


    /*Saves estimated signal, true signal and noisy signal on a txt file*/
    ofstream outfile1 ("signals.txt");
    outfile1 << "estimated_signal" << " " << "true_signal" << " " << "noisy_signal" << endl;
    for (nat i = 0; i < cols(states);i++){
      outfile1 << get_data(estimated_signal,0,i) << " " 
      << get_data(signal,0,i) << " " 
      <<  get_data(y,0,i) << endl;
    }
    outfile1.close();

    /*Saves states evolution on a txt file*/
    ofstream outfile2 ("states.txt");
    outfile2 << "A" << " " << "Xi"  << " " << "w" << endl;
    for (nat i = 0; i < cols(states);i++){
      outfile2 << get_data(states,0,i) << " " 
      << get_data(states,1,i) << " "
      << get_data(states,2,i) << endl;
    }
    outfile2.close();

    
    /*Delete all matrices to avoid memory leaks*/
    delete_matrix(out_matrix);
    delete_matrix(states);
    delete_matrix(estimated_signal);
    delete_matrix(I);
    delete_matrix(P_prior);
    delete_matrix(Q);
    delete_matrix(R);
    delete_matrix(y);
    delete_matrix(x_prior);
    delete_matrix(signal);
    delete_EKF(yk);
    delete_EKF(xk);
    delete_EKF(G);
    delete_EKF(C);
    delete_EKF(Phi);
    delete_EKF(U);

}



