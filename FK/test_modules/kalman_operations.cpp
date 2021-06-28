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
  
  nat n_samples = 150; 
  nat n = 2, m = 1;
  TMatrix P = create_matrix(n,n);
  TMatrix R = create_matrix(m,m);
  TMatrix C = create_matrix(m,n);
  TMatrix xk = create_matrix(n,1);
  TMatrix y = create_matrix(m,n_samples);
  TMatrix signal = create_matrix(m,n_samples);
  TMatrix Q = create_matrix(n,n);
  TMatrix Phi = create_matrix(n,n);
  TMatrix I = create_matrix(n,n);

  insert_data(P,0,0,10);
  insert_data(P,1,0,0);
  insert_data(P,0,1,0);
  insert_data(P,1,1,10);
  insert_data(R,0,0,9.0); 
  insert_data(C,0,0,1);
  insert_data(C,0,1,0);
  insert_data(xk,0,0,0.1);
  insert_data(xk,1,0,0.5);

  insert_data(Q,0,0,0);
  insert_data(Q,1,0,0);
  insert_data(Q,0,1,0);
  insert_data(Q,1,1,0);

  insert_data(Phi,0,0,1.0);
  insert_data(Phi,0,1,1.0);
  insert_data(Phi,1,0,0.0);
  insert_data(Phi,1,1,1.0);
  
 
  insert_data(I,0,0,1);
  insert_data(I,1,0,0);
  insert_data(I,0,1,0);
  insert_data(I,1,1,1); 


  const double mean = 0.0;
  const double stddev = 3;
  std::default_random_engine generator;
  std::normal_distribution<double> dist(mean, stddev);
  
  for (nat i = 0; i < cols(y);i++){ 
    insert_data(signal,0,i,0.065*i);
    insert_data(y,0,i,0.065*i+dist(generator));
  }

  auto start = chrono::steady_clock::now();
  
  
  TMatrix out_matrix = Kalman_filter(I,P,Q,Phi,xk,y,C,R);

  auto end = chrono::steady_clock::now();
   cout << "Elapsed time in nanoseconds: " 
	<< chrono::duration_cast<chrono::nanoseconds>(end - start).count()
        << " ns" << endl;
   
  TMatrix p_signal = get_signal(out_matrix, rows(y));
  TMatrix states = get_state(out_matrix,rows(y));




  ofstream outfile1 ("signals.txt");
  outfile1 << "p_signal" << " " << "signal" << " " << "noisy_signal" << endl;
  for (nat i = 0; i < cols(states);i++){
	    outfile1 << get_data(p_signal,0,i) << " " 
	    << get_data(signal,0,i) << " " 
	    <<  get_data(y,0,i) << endl;
  }
  outfile1.close();


  ofstream outfile2 ("states.txt");
  outfile2 << "rk" << " " << "pk" << endl;
  for (nat i = 0; i < cols(states);i++){
    outfile2 << get_data(states,0,i) << " " 
    << get_data(states,1,i) << endl;
  }
  outfile2.close();
  


  delete_matrix(p_signal);
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



