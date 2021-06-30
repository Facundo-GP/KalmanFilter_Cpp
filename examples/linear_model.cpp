#include <stdio.h>
#include "../include/matrix.h"
#include "../include/kalman.h"
#include <math.h>
#include <iostream>
#include <chrono>
#include <unistd.h>
#include <random>
#include <fstream>

using namespace std;

int main()

{
  
  /*n_samples = number of iterations*/
  nat n_samples = 150; 
  
  /*n = states size, m=output size*/
  nat n = 2, m = 1;
  
  /*P*/
  TMatrix P = create_matrix(n,n);
  double P_data[4] = {1000,0,0,1000};
  fill_matrix(P,P_data);
  
  /*R*/
  TMatrix R = create_matrix(m,m);
  double R_data[1] = {9.0};
  fill_matrix(R,R_data);

  /*C*/
  TMatrix C = create_matrix(m,n);
  double C_data[2] = {1.0,0};
  fill_matrix(C,C_data);

  /*xk prior*/
  TMatrix xk = create_matrix(n,1);
  double xk_data[2] = {0.1,0.5};
  fill_matrix(xk,xk_data);
  
  /*Q, default zero matrix*/
  TMatrix Q = create_matrix(n,n);

  /*Phi*/
  TMatrix Phi = create_matrix(n,n);
  double Phi_data[4] = {1.0,1.0,0,1.0};
  fill_matrix(Phi,Phi_data);

  /*I*/
  TMatrix I = identity(n); 
  
  /*Matrix that will contain true signal for comparation*/
  TMatrix signal = create_matrix(m,n_samples);
  
  /*y*/
  TMatrix y = create_matrix(m,n_samples);
 
  /*noise*/
  const double mean = 0.0;
  const double stddev = 3;
  std::default_random_engine generator;
  std::normal_distribution<double> dist(mean, stddev);
  
  /*Signal generation*/
  for (nat i = 0; i < cols(y);i++){ 
    insert_data(signal,0,i,0.065*i);
    insert_data(y,0,i,0.065*i+dist(generator));
  }

  
  /*Kalman filter*/
  TMatrix out_matrix = Kalman_filter(I,P,Q,Phi,xk,y,C,R);

  /*States and predicted signal evolution*/
  TMatrix estimated_signal = get_signal(out_matrix, rows(y));
  TMatrix states = get_state(out_matrix,rows(y));


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
  outfile2 << "rk" << " " << "pk" << endl;
  for (nat i = 0; i < cols(states);i++){
    outfile2 << get_data(states,0,i) << " " 
    << get_data(states,1,i) << endl;
  }
  outfile2.close();
  


  /*Delete all matrices to avoid memory leaks*/
  delete_matrix(estimated_signal);
  delete_matrix(signal);
  delete_matrix(states);
  delete_matrix(out_matrix);
  delete_matrix(P);
  delete_matrix(R);
  delete_matrix(C);
  delete_matrix(xk);
  delete_matrix(y);
  delete_matrix(Q);
  delete_matrix(Phi);
  delete_matrix(I);


   
}



