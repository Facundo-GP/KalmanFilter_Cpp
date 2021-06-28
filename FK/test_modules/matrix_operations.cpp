#include <stdio.h>
#include <assert.h>
#include "include/matrix.h"
#include <math.h>

int main(){
  int n = 5, m = 5;
  TMatrix A = create_matrix(n,m);
  TMatrix B = create_matrix(m,m);
  for (int i = 0; i < n;i++){
    for (int j = 0; j < m; j++){
      insert_data(A,i,j,static_cast <double> (rand()) / static_cast <double> (RAND_MAX/10));
    }
  }
  /*
  double datos[3][3] ;
  datos[0][0] = 0.5;
  datos[0][1] = -1.3;
  datos[0][2] = 3;
  datos[1][0] = 2;
  datos[1][1] = -10; 
  datos[1][2] = 15;
  datos[2][0] = -12.3;
  datos[2][1] = 0;
  datos[2][2] = 2;
  for (int i = 0; i < m;i++){
    for (int j = 0; j < m; j++){
      insert_data(B,i,j,datos[i][j]);
    }
  }*//*
  printf("Matriz A: \n");
  display_matrix(A);*/
  /*printf("\n Matriz B: \n ");
  display_matrix(B);*/

 /* printf("\n Matriz A + B : \n");
  TMatrix suma = sum_matrices(A,B);
  display_matrix(suma);
  delete_matrix(suma)*/
  /*
  printf("\n AxB : \n");
  TMatrix mult = multiply(A,B);
  display_matrix(mult);
  delete_matrix(mult);
  
  printf("\n AxB transpuesta: \n");
  TMatrix trans = transpose(B);
  display_matrix(trans);
  delete_matrix(trans);


  printf("\n (A.TxA + B).T + B \n");
  TMatrix AT = transpose(A);
  TMatrix ATA = multiply(AT,A);
  delete_matrix(AT);
  TMatrix ATAPB = sum_matrices(ATA,B);
  delete_matrix(ATA);
  TMatrix ATAPBT = transpose(ATAPB);
  delete_matrix(ATAPB);
  TMatrix res = sum_matrices(ATAPBT,B);
  display_matrix(res);
 */
  /*printf("\n Determinante: \n");
  double determinante = det(A);
  printf("%f \n",determinante);*/
/*
  printf("\nSub matriz\n");
  TMatrix sub = submatrix(B,1,1);
  display_matrix(sub);

  
  for (nat k = 0; k < 1000;k++){
    printf("\n Inversa \n");
    TMatrix inv = inverse(A);
    delete_matrix(inv);
  }*/
  
  printf("Sub_column: \n");
  TMatrix S = sub_column(A,4);
  display_matrix(A);
  printf("\n\n");
  display_matrix(S);
  fill_column(A,0,S);
  display_matrix(A);
  /*display_matrix(inv);*/

  /*delete_matrix(inv);*/
  delete_matrix(B);
  delete_matrix(A);
  /*delete_matrix(sub);
  delete_matrix(res);
  delete_matrix(ATAPBT);*/


}
