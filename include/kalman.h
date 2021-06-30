#include "matrix.h"

typedef struct _rep_kalman_f* TEKF_func;


/*Returns a TEKF_func type used to describe C, F, G,U ,xk and yk evolution over iterations of EKF*/
/*len is the number of states*/
/*Precondition: params has len elements*/
TEKF_func create_kalman_f(TMatrix params, nat iter_param,TMatrix (*f)(TMatrix,nat));


void replace_f_params(TEKF_func EKF_f , TMatrix params, nat iter_params);


TMatrix evaluate_function(TEKF_func EKF_f);


void delete_EKF(TEKF_func EKF_f);



TMatrix get_params(TEKF_func EKF_f);



TMatrix EKF_gain(TMatrix P, TMatrix C, TMatrix R, TMatrix U);


TMatrix EKF_update_state(TMatrix xk, TMatrix K, TMatrix yk,TMatrix y_pred);

TMatrix EKF_P_projection(TMatrix P, TMatrix Phi, TMatrix G,TMatrix Q);


/*Return the Kalman gain K*/
/*Preconditions: P nxn, C mxn, R mxm*/
TMatrix Kalman_gain(TMatrix P, TMatrix C, TMatrix R);

/*Update state estiamtion*/
/*Preconditions: xk nx1, K nxm, yk mx1 C mxn*/
TMatrix update_state(TMatrix xk, TMatrix K, TMatrix yk, TMatrix C);


/*Update P estimation*/
/*Preconditions I nxn , P nxn, K nxm,  C mxn */
TMatrix update_P(TMatrix I, TMatrix P, TMatrix K, TMatrix C);


/*State projection*/
/*Preconditions: x_updated nx1, Phi nxn*/
TMatrix state_projection(TMatrix x_updated, TMatrix Phi);


/*P projection*/
/*Precondition: P nxn, Q nxn, Phi nxn*/
TMatrix P_projection(TMatrix P, TMatrix Q, TMatrix Phi);


/*Return the kalman filter output for yk samples */


/*Precondition: I nxn, P nxn, Q nxn, Phi nxn, xk nx1 , yk  m x num_samples, R mxm 
              C mxn num_samples > 0*/


TMatrix Kalman_filter(TMatrix I, TMatrix P, TMatrix Q, TMatrix Phi, 
		     TMatrix xk, TMatrix y, TMatrix C,
		     TMatrix R);



/*Returns the signal estimated evolution of the Kalman_filter function output matrix*/
/*Preconditions:  0 < signal_out_dim < rows(KF_out)-1*/
TMatrix get_signal(TMatrix KF_out, nat signal_rows_dim);

/*Returns the state estimated evolution of the Kalman_filter function output matrix*/
/*Preconditions: 0 < signal_out_dim <  rows(KF_out)-1*/
TMatrix get_state(TMatrix KF_out, nat signal_rows_dim);



/*Return the EKF output for yk samples */

/*C and xk are functions that implement the update rule of the non linear problem i.e the first order aproximation used to aproximate C and xk*/
/*double* params are the params of the function and nat the len of double* */

/*Precondition: I nxn, P nxn, Q nxn, Phi nxn, xk nx1 , yk  m x num_samples, R mxm 
              C mxn num_samples > 0*/


TMatrix EKF(TMatrix I, TMatrix P_prior, TMatrix Q, TMatrix R,TMatrix y,
	    TMatrix x_prior, TEKF_func C,
	    TEKF_func Phi, TEKF_func G, TEKF_func U, TEKF_func xk,
	    TEKF_func yk);
