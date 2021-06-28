

typedef unsigned int nat;

/*Module that declarate matrix operations*/

typedef struct _rep_matrix * TMatrix;


/*Create an empty matrix of size with n rows and m columns
Precondition: n, m > 0 O(n) complexity  */
TMatrix create_matrix(nat n, nat m);


/*Delete the matrix and all its components O(1) complexity*/
void delete_matrix(TMatrix A);


/*Num of rows of a matrix*/
nat rows(TMatrix A);

/*Num of cols of a matrix*/
nat cols(TMatrix A);

/*Insert a number in the i rown and j column  O(1) complexity */
/*Precondition: i and j are in the allowed matrix index*/
void insert_data(TMatrix A,nat row, nat col,double data);

/*Prints a matrix O(n*m) complexity*/
void display_matrix(TMatrix A);


/*Returns  the transpose of matrix O(n*m) complexity*/
/*TMatrix transpose(TMatrix matrix);*/
TMatrix transpose(TMatrix A);


/*Returns the matrix sum between A matrix and B matrix*/
/*Precondition: shapes of matrices are equal*/
TMatrix sum_matrices(TMatrix A, TMatrix B);



/*Return subtraction betweehn A and B (A-B)*/
/*Preconditions: shapes of matrices are equal*/
TMatrix substract_matrices(TMatrix A, TMatrix B);


/*Returns the matrix multiplication between A matrix and B matrix*/
/*Precondition: A columns number and B rows number are equal*/
TMatrix multiply(TMatrix A, TMatrix B);


/*Swap i and j rows if row == true or columns if rows == false of matrix A with size nxm */
/*Precondition: if i and j are rows, j <= n-1, i <= n-1*/
/*if i and j are columns, j <= m-1, i <= m-1*/
void swap_rows_cols(TMatrix A, nat i, nat j,bool rows);


/*Copy data of matrix A in another*/
TMatrix copy(TMatrix A);

/*Returns  i column of matrix A*/
/*Precondition: 0 < i < cols(A)*/
TMatrix sub_column(TMatrix A, nat i);

/*Fills i-th column of A with data on D matrix */
/*Precondition, D is a matrix with  one column and rows(A) rows*/
void fill_column(TMatrix A, nat i, TMatrix D);

/*Returns the submatrix of A withouth i-th column and row*/
/*Precondition: i a valid index of matrix A, square matrix*/
/*A shape greather than 1x1*/
TMatrix submatrix(TMatrix A, nat row_index,nat col_index);

/*Retuns the determinant of Matrix A*/
/*Precondition: non singular matrix, square matrix*/
double det(TMatrix A);

/*Return inverse of the matrix*/
/*Precondition: Non singular matrix, square matrix*/
TMatrix inverse(TMatrix A);

/*Fills the sub matrix of A delimited by r1 and r2 rows & by c1 & c2 columns with D matrix data */
/*Precondition: r1 <= r2, c1 <= c2  r1,c1 >= 0  r1 < rows(A), c1 < cols(A) 
rows(D) = r2-r1, cols(D) = c2-c1*/
void fill_sub(TMatrix A, nat r1, nat r2, nat c1, nat c2,TMatrix D);


/*Returns data on i-th row and j-th columns of matrix A*/
/*Preconditions, 0 <= i < rows(A)   0 <= j < cols(A)*/
double get_data(TMatrix A, nat i, nat j);



/*Returns the sub matrix of A delimited by r1 and r2 rows & by c1 & c2 columns with D matrix data */
/*Precondition: r1 <= r2, c1 <= c2  r1,c1 >= 0  r1 < rows(A), c1 < cols(A)*/
TMatrix sub_matrix(TMatrix A, nat r1, nat r2, nat c1, nat c2);

/*Fills matrix A by row with new_data */
/**Preconditions: new_data has exactly rows(A)*cols(A) elems*/
void fill_matrix(TMatrix A, double* new_data);


/*Return the nxn identity*/
TMatrix identity(nat n);
