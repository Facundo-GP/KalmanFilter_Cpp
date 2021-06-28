#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "../include/matrix.h"


struct _rep_matrix{
  nat rows,columns;
  double ** data;

};


TMatrix create_matrix(nat n, nat m)

{
    TMatrix matrix  = new _rep_matrix;
    matrix->rows = n;
    matrix->columns = m;
    matrix->data = new double*[n];
    for (nat i = 0;i < n;i++){
      matrix->data[i] = new double[m];
      for (nat j = 0; j < m; j++){
        matrix->data[i][j] = 0.0;
      }
    }
    
    return matrix;

}

void delete_matrix(TMatrix A)

{
    for (nat i = 0; i < A->rows;i++){
      delete [] A->data[i];
    }
    delete [] A->data;
    delete A;
}

nat rows(TMatrix A)

{
    return A->rows;
}

nat cols(TMatrix A)

{
    return A->columns;
}


void insert_data(TMatrix A,nat row, nat col,double data)

{
    assert(row <= A->rows && col <= A->columns); 
    A->data[row][col] = data;
  
}


void display_matrix(TMatrix A)

{

    for (nat i = 0;i < A->rows;i++){
      for (nat j = 0; j < A->columns;j++){
        if (j == 0 && A->columns != 1){
          printf("| %.4f ,",A->data[i][j]);
        }
	else if ( j == A->columns-1 && A->columns != 1){
          printf(" %.4f |",A->data[i][j]);
        }
	else if (A->columns == 1){
          printf("| %.4f  |",A->data[i][j]);
	}
        else{printf(" %.4f ,",A->data[i][j]);}
      }
      printf("\n");
    }

}


TMatrix transpose(TMatrix A)

{
     TMatrix transpose_matrix = create_matrix(A->columns,A->rows);
      
    /*Case of square matrix*/
    if (A->rows == A->columns){
      double aux;
      for (nat i = 0; i < A->rows;i++){
        for (nat j = 0; j < A->columns;j++){
          aux = A->data[j][i];
          transpose_matrix->data[j][i] = A->data[i][j];
	  transpose_matrix->data[i][j] = aux;
        }
      }
    }
    
    
    /*Non square matrix*/
    else {
        for (nat i =0; i < A->rows;i++){
          for (nat j = 0; j < A->columns;j++){
	    insert_data(transpose_matrix,j,i,A->data[i][j]);
	  }
        }
    }
    return transpose_matrix;

}


TMatrix sum_matrices(TMatrix A, TMatrix B)
{

    assert(A->rows == B->rows && A->columns == B->columns);
    TMatrix C = create_matrix(A->rows,A->columns);

    for (nat i =0; i < C->rows;i++){
      for (nat j = 0; j < C->columns;j++){
	  insert_data(C,i,j,A->data[i][j]+B->data[i][j]);
       }
    }
    return C;
}

TMatrix substract_matrices(TMatrix A, TMatrix B)


{

    assert(A->rows == B->rows && A->columns == B->columns);
    TMatrix neg_B = copy(B);
    for (nat i =0; i < neg_B->rows;i++){
      for (nat j = 0; j < neg_B->columns;j++){
	  insert_data(neg_B,i,j,-B->data[i][j]);
       }
    }

    TMatrix C = sum_matrices(A,neg_B);
    delete_matrix(neg_B);
    return C;
    
} 


void swap_rows_cols(TMatrix A,nat i,nat j,bool rows)


{ 

  if (rows){
    assert((i < A->rows) && (j < A->rows));
    double aux;
    for (nat k = 0; k < A->columns;k++){
      aux = A->data[i][k];
      A->data[i][k] = A->data[j][k];
      A->data[j][k] = aux;
    }
  }

  else{
    assert((i < A->columns) && (j < A->columns));
    double aux;
    for (nat k = 0; k < A->rows;k++){
      aux = A->data[k][i];
      A->data[k][i] = A->data[k][j];
      A->data[k][j] = aux;
    }
  }
}

TMatrix multiply(TMatrix A,TMatrix B)

{
  assert(A->columns == B->rows);
  TMatrix C = create_matrix(A->rows,B->columns);
  for (nat i = 0; i < A->rows; i++){
    for (nat j = 0; j < B->columns; j++){
      for (nat k = 0; k < A->columns;k++){
        C->data[i][j] += A->data[i][k]*B->data[k][j];
      }
    }
  }
  return C;

}


TMatrix copy(TMatrix A)

{
  TMatrix A_copy = create_matrix(A->rows,A->columns);
  for (nat i = 0; i < A->rows;i++){
    for (nat j = 0; j < A->columns; j++){
      A_copy->data[i][j] = A->data[i][j];
    }
  }
  return A_copy;
}

double det(TMatrix A)

{

    assert(A->rows == A->columns);/*square matrix*/

    TMatrix A_copy = copy(A);
    double d = 1;
    
    if (A->rows == 1){
      d = A_copy->data[0][0];
      delete_matrix(A_copy);
      return d;
    }
    else if (A->rows == 2){
      d = A_copy->data[0][0]*A_copy->data[1][1]-A_copy->data[1][0]*A_copy->data[0][1];
      delete_matrix(A_copy);
      return d;
    }

    else {
      for (nat i = 0; i < A_copy->rows;i++){
        if(A_copy->data[i][i] == 0.0){
	  
	  if (i == A_copy->rows-1){ /*singular matrix*/
	    delete_matrix(A_copy);
	    return 0.0;
	  }
	  else{
            swap_rows_cols(A_copy,i,i+1,true);
	  }
	}
	for (nat j = i+1; j < A_copy->rows;j++){
          double ratio = A_copy->data[j][i]/A_copy->data[i][i];
	  for (nat k = 0; k < A_copy->rows;k++){
	    A_copy->data[j][k] -= ratio*A_copy->data[i][k];
	  }
	}
      }
      for (nat i = 0; i < A_copy->rows;i++){
        d *= A_copy->data[i][i];
      }
    }

    delete_matrix(A_copy); 
    return d;

}


TMatrix sub_column(TMatrix A, nat i)

{
   assert(i < cols(A));
   assert(i >= 0);

   TMatrix row_matrix = create_matrix(rows(A),1);

   for (nat j = 0; j < rows(A);j++){
     insert_data(row_matrix,j,0,A->data[j][i]);
   }

   return row_matrix;
}


void fill_column(TMatrix A, nat i, TMatrix D)

{
    assert(rows(A) == rows(D));
    assert(cols(D) == 1);

   for (nat j = 0; j < rows(A);j++){
     insert_data(A,j,i,D->data[j][0]);
   }

}

TMatrix submatrix(TMatrix A,nat row_index,nat col_index)

{
    assert(A->rows == A->columns && A->rows > 1);
    TMatrix sub = create_matrix(A->rows-1,A->columns-1);
    for (nat i = 0; i < A->rows;i++){
      for (nat j = 0; j < A->columns; j++){
        if (i <  row_index && j < col_index){
          sub->data[i][j] = A->data[i][j];
        } 
        else if (i >  row_index && j < col_index){
          sub->data[i-1][j] = A->data[i][j];
        }
        else if (i <  row_index && j > col_index){
          sub->data[i][j-1] = A->data[i][j];
        }
        else if (i >  row_index && j > col_index){
    	  sub->data[i-1][j-1] = A->data[i][j];
      
        }

      }
    }
    return sub;

}


TMatrix inverse(TMatrix A)

{

  TMatrix inv = copy(A);
  double d  = det(A);
  assert(d != 0.0);
  if (A->rows == 1){
    inv->data[0][0] = 1/d;
    return inv;
  }
  else{
    for (nat i = 0; i < A->rows; i++){
      for (nat j = 0; j < A->columns; j++){
        TMatrix sub = submatrix(A,i,j);
        inv->data[i][j] = pow(-1.0,i+j)*det(sub);
        delete_matrix(sub);
      }
    }
  }

  TMatrix inv_t = transpose(inv);
  delete_matrix(inv);
   
  for (nat i = 0; i < A->rows; i++){
    for (nat j = 0; j < A->columns; j++){
      inv_t->data[i][j] /= det(A);
    }
  }
  

  return inv_t;

}

void fill_sub(TMatrix A, nat r1, nat r2, nat c1, nat c2, TMatrix D)

{
    assert( (r1 >= 0) && (c1 >= 0));
    assert((r2 <= rows(A)) && (c2 <= cols(A)) );
    assert((r2 >= r1) && (c2 >= c1));
    assert((r2-r1 == rows(D)) && (c2-c1 == cols(D)));


    for (nat i = r1; i <= r2;i++){
      for (nat j = c1; j <= c2; j++){
        A->data[i][j] == D->data[i-r1][j-c1];
      }
    }
}



double get_data(TMatrix A, nat i, nat j)

{
    assert((0 <= i) && (i < rows(A)));
    assert((0 <= j ) && (j < cols(A)));
    return A->data[i][j];

}

void fill_matrix(TMatrix A, double* new_data)

{
  
  nat cont = 0;
  for (nat i = 0; i < rows(A);i++){
    for (nat j = 0; j < cols(A);j++){
      cont++;
      A->data[i][j] = new_data[cont-1];
    }
  }
  

}

TMatrix sub_matrix(TMatrix A, nat r1, nat r2, nat c1, nat c2)

{
    assert( (r1 >= 0) && (c1 >= 0));
    assert((r2 <= rows(A)) && (c2 <= cols(A)) );
    assert((r2 >= r1) && (c2 >= c1));
    
    
    TMatrix sub = create_matrix(r2-r1+1,c2-c1+1);
    for (nat i = r1; i <= r2;i++){
      for (nat j = c1; j <= c2; j++){
        sub->data[i-r1][j-c1] = A->data[i][j];
      }
    }
    return sub;
}

TMatrix identity(nat n)

{
  TMatrix I = create_matrix(n,n);
  for (nat i = 0; i < n; i++){
    insert_data(I,i,i,1);
  }
  return I;

}

