#include <stdio.h>
#include <assert.h>
#include "../include/kalman.h"
#include "../include/matrix.h"
#include <cstring>



struct _rep_kalman_f{

   TMatrix (*function)(TMatrix ,nat);
   TMatrix state_params;
   nat iter_param;

};

TEKF_func create_kalman_f(TMatrix params, nat iter_param,TMatrix (*f)(TMatrix,nat))

{
  TEKF_func EKF_f = new _rep_kalman_f;
  EKF_f->function = f;
  EKF_f->state_params = create_matrix(rows(params),1);
  
  for (nat i = 0; i < rows(params);i++){
    insert_data(EKF_f->state_params,i,0,get_data(params,i,0));
  }

  EKF_f->iter_param = iter_param;

  return EKF_f;

}


void replace_f_params(TEKF_func EKF_f, TMatrix params, nat iter_param)

{
  assert(rows(params) == rows(EKF_f->state_params));
  for (nat i = 0; i < rows(params);i++){
    insert_data(EKF_f->state_params,i,0,get_data(params,i,0));
  }
  EKF_f->iter_param = iter_param;

}

TMatrix evaluate_function(TEKF_func EKF_f)


{

  TMatrix res = EKF_f->function(EKF_f->state_params,EKF_f->iter_param);
  return res;

}

void delete_EKF(TEKF_func EKF_f)

{
  delete_matrix(EKF_f->state_params);
  delete EKF_f;

}



TMatrix get_params(TEKF_func EKF_f)

{
  return EKF_f->state_params;
}



TMatrix EKF_gain(TMatrix P, TMatrix C, TMatrix R,TMatrix U)

{
    /*Kalman gain*/
    TMatrix C_T = transpose(C);
    TMatrix PxC_T = multiply(P,C_T);
    TMatrix CxP = multiply(C,P);
    TMatrix CxPxC_T = multiply(CxP,C_T);

    TMatrix UxR = multiply(U,R);
    TMatrix U_T = transpose(U);
    TMatrix UxRxU_T = multiply(UxR,U_T);
    TMatrix CxPxC_TpUxRxU_T = sum_matrices(CxPxC_T,UxRxU_T);
    TMatrix CxPxC_TpUxRxU_T_I = inverse(CxPxC_TpUxRxU_T);
    TMatrix K = multiply(PxC_T,CxPxC_TpUxRxU_T_I);
    
   

    /*Delete intermediate matrices*/
    delete_matrix(C_T);
    delete_matrix(PxC_T);
    delete_matrix(CxP);
    delete_matrix(CxPxC_T);
    delete_matrix(UxR);
    delete_matrix(U_T);
    delete_matrix(UxRxU_T);
    delete_matrix(CxPxC_TpUxRxU_T);
    delete_matrix(CxPxC_TpUxRxU_T_I);
    
    return K;



}



TMatrix EKF_update_state(TMatrix xk, TMatrix K, TMatrix yk, TMatrix y_pred)

{
 
    /*update state vector*/
    TMatrix ykmy_pred = substract_matrices(yk,y_pred);
    TMatrix Kykmy_pred = multiply(K,ykmy_pred);
    TMatrix xk_updated = sum_matrices(xk,Kykmy_pred);

    /*delete intermediate matrices*/
    delete_matrix(ykmy_pred);
    delete_matrix(Kykmy_pred);

    return xk_updated;

 
}


TMatrix EKF_P_projection(TMatrix P, TMatrix Phi, TMatrix G, TMatrix Q)

{

    /*P projection*/
    TMatrix Phi_T = transpose(Phi);
    TMatrix PhixP = multiply(Phi,P);
    TMatrix PhixPxPhi_T = multiply(PhixP,Phi_T);

    TMatrix G_T = transpose(G);
    TMatrix GxQ = multiply(G,Q);
    TMatrix GxQxG_T = multiply(GxQ,G_T);

    TMatrix P_projected = sum_matrices(PhixPxPhi_T,GxQxG_T);
    
    /*delete used matrices*/
    delete_matrix(Phi_T);
    delete_matrix(PhixP);
    delete_matrix(PhixPxPhi_T);
    delete_matrix(G_T);
    delete_matrix(GxQ);
    delete_matrix(GxQxG_T);

    return P_projected;


}


TMatrix Kalman_gain(TMatrix P, TMatrix C, TMatrix R)

{
    /*Preconditions*/
    assert(rows(P) == cols(P)); 
    assert(rows(R) == cols(R)); 
    assert(cols(C) == rows(P));
    assert(rows(C) == rows(R));
 
    /*Kalman gain*/
    TMatrix C_T = transpose(C);
    TMatrix PxC_T = multiply(P,C_T);
    TMatrix CxP = multiply(C,P);
    TMatrix CxPxC_T = multiply(CxP,C_T);
    TMatrix CxPxC_TpR = sum_matrices(CxPxC_T,R);
    TMatrix CxPxC_TpR_I = inverse(CxPxC_TpR);
    TMatrix K = multiply(PxC_T,CxPxC_TpR_I);
    
     

    /*Delete intermediate matrices*/
    delete_matrix(C_T);
    delete_matrix(PxC_T);
    delete_matrix(CxP);
    delete_matrix(CxPxC_T);
    delete_matrix(CxPxC_TpR);
    delete_matrix(CxPxC_TpR_I);
    
    return K;

}



TMatrix update_state(TMatrix xk, TMatrix K, TMatrix yk, TMatrix C)

{

    /*Preconditions*/
    assert(cols(xk) == 1);
    assert(cols(yk) == 1);
    assert(cols(K) == rows(C));
    assert(rows(K) == cols(C));
    assert(cols(K) == rows(yk));
    assert(rows(K) == rows(xk));
    
    /*update state vector*/
    TMatrix Cxxk = multiply(C,xk);
    TMatrix ykmCxxk = substract_matrices(yk,Cxxk);
    TMatrix KykmCxxk = multiply(K,ykmCxxk);
    TMatrix xk_updated = sum_matrices(xk,KykmCxxk);

    /*delete intermediate matrices*/
    delete_matrix(Cxxk);
    delete_matrix(ykmCxxk);
    delete_matrix(KykmCxxk);

    return xk_updated;

}

TMatrix update_P(TMatrix I,TMatrix P,TMatrix K,TMatrix C)

{
    
    /*Preconditions*/
    assert(rows(P) == cols(P));
    assert(cols(K) == rows(C));
    assert(rows(K) == cols(C));
    assert(rows(P) == rows(K));
    assert( (rows(I) == cols(I)) && (rows(I) == rows(K)));

    /*update P*/
    TMatrix KxC = multiply(K,C);
    TMatrix ImKxC = substract_matrices(I,KxC);
    TMatrix P_updated = multiply(ImKxC,P);

    /**delete used matrices*/
    delete_matrix(KxC);
    delete_matrix(ImKxC);

    return P_updated;

}


TMatrix state_projection(TMatrix x_updated,TMatrix Phi)

{


    /*Preconditions*/
    assert(cols(x_updated) == 1);
    assert(rows(x_updated) == rows(Phi));
    assert(rows(Phi) == cols(Phi));

    /*state projection*/
    TMatrix x_projected = multiply(Phi,x_updated);
    return x_projected;

}


TMatrix P_projection(TMatrix P, TMatrix Q, TMatrix Phi)

{
    
    /*Preconditions*/
    assert(rows(P) == cols(P));
    assert(rows(Q) == cols(Q));
    assert(rows(Phi) == cols(Phi));
    assert(rows(Q) == rows(Phi));
    assert(rows(Q) == rows(P));
    /*P projection*/
    TMatrix Phi_T = transpose(Phi);
    TMatrix PhixP = multiply(Phi,P);
    TMatrix PhixPxPhi_T = multiply(PhixP,Phi_T);
    TMatrix P_projected = sum_matrices(PhixPxPhi_T,Q);
    
    /*delete used matrices*/
    delete_matrix(Phi_T);
    delete_matrix(PhixP);
    delete_matrix(PhixPxPhi_T);

    return P_projected;

}


TMatrix Kalman_filter(TMatrix I,TMatrix P, TMatrix Q, TMatrix Phi, TMatrix xk,
		      TMatrix y, TMatrix C, TMatrix R)

{
    assert(rows(P) == cols(P));
    assert(rows(Q) == cols(Q));
    assert(rows(Phi) == rows(Phi));
    assert(rows(R) == cols(R));
    assert((rows(I) == cols(I)) && (rows(I) == rows(P)));
    assert( (rows(P) == rows(Q)) && (rows(P) == rows(Phi)) && (rows(R) == rows(C)));
    assert( (rows(xk) == cols(C)) && (rows(y) == rows(C) ));
    assert( (cols(xk) == 1) ); 
 

    TMatrix xk_c = copy(xk);
    TMatrix P_c = copy(P);
   
    
      
    TMatrix out_matrix = create_matrix(rows(xk)+rows(y),cols(y));
      
    TMatrix K = Kalman_gain(P_c,C,R);
    TMatrix yk = sub_column(y,0);
    TMatrix xk_updated = update_state(xk_c,K,yk,C);
    TMatrix P_updated = update_P(I,P_c,K,C);
    TMatrix x_projected = state_projection(xk_updated,Phi);
    TMatrix P_projected = P_projection(P_updated,Q,Phi);
    TMatrix y_pred = multiply(C,xk_updated);
	
    for(nat j = 0; j < rows(out_matrix);j++){ 
      if (j < rows(y)){
        insert_data(out_matrix,j,0,get_data(y_pred,j,0));
      } 
      else{
        insert_data(out_matrix,j,0,get_data(xk_updated,j-rows(y),0));
      }
    }

    delete_matrix(y_pred);
    delete_matrix(yk);
    delete_matrix(K);
    delete_matrix(xk_updated);
    delete_matrix(P_updated);  
    delete_matrix(P_c);
    delete_matrix(xk_c);
        
    P_c = P_projected;
    xk_c = x_projected;


    for (nat i = 1; i < cols(y); i++) {
      K = Kalman_gain(P_c,C,R);
      yk = sub_column(y,i);
      xk_updated = update_state(xk_c,K,yk,C);
      P_updated = update_P(I,P_c,K,C);  
      x_projected = state_projection(xk_updated,Phi);
      P_projected = P_projection(P_updated,Q,Phi); 
      y_pred = multiply(C,xk_updated);
	
      for(nat j = 0; j < rows(out_matrix);j++){ 
        if (j < rows(y)){
          insert_data(out_matrix,j,i,get_data(y_pred,j,0));
	}
	else{
          insert_data(out_matrix,j,i,get_data(xk_updated,j-rows(y),0));
        }
      }

      delete_matrix(y_pred);
      delete_matrix(yk);
      delete_matrix(K);
      delete_matrix(xk_updated);
      delete_matrix(P_updated);
      delete_matrix(P_c);
      delete_matrix(xk_c);
        
      P_c = P_projected;
      xk_c = x_projected;

    }
    delete_matrix(P_c);
    delete_matrix(xk_c);
    return out_matrix;

}



TMatrix get_signal(TMatrix KF_out, nat signal_out_dim)

{ 
  assert( (0 < signal_out_dim) && (signal_out_dim < rows(KF_out)-1));
  TMatrix signal = sub_matrix(KF_out,0,signal_out_dim-1,0,cols(KF_out)-1);
  return signal;

}

TMatrix get_state(TMatrix KF_out, nat signal_out_dim)

{ 
  assert( (0 < signal_out_dim) && (signal_out_dim < rows(KF_out)-1));
  TMatrix signal = sub_matrix(KF_out,signal_out_dim,rows(KF_out)-1,0,cols(KF_out)-1);
  return signal;

}

TMatrix EKF(TMatrix I, TMatrix P_prior, TMatrix Q, TMatrix R,TMatrix y,TMatrix x_prior,
            TEKF_func C,TEKF_func Phi,TEKF_func G, 
	   TEKF_func U, TEKF_func xk, TEKF_func yk)

{


    TMatrix xk_c = copy(x_prior);
    TMatrix P_c = copy(P_prior);
    TMatrix out_matrix = create_matrix(rows(xk_c)+rows(y),cols(y));
    

    TMatrix C_k = evaluate_function(C);
    TMatrix U_k = evaluate_function(U);
    /*
    printf("C: \n");
    display_matrix(C_k);
     
    printf("P: \n");
    display_matrix(P_c);
*/
    
    TMatrix K = EKF_gain(P_c,C_k,R,U_k);
    TMatrix y_pred = evaluate_function(yk);
    TMatrix yk_m = sub_column(y,0);
/*
    printf("K: \n");
    display_matrix(K);
  */  TMatrix xk_updated = EKF_update_state(xk_c,K,yk_m,y_pred); 
    TMatrix P_updated = update_P(I,P_c,K,C_k);
    

    /*projections*/

    replace_f_params(xk,xk_updated,0);
    TMatrix x_projected = evaluate_function(xk);
    
    replace_f_params(G,xk_updated,0);
    TMatrix G_k = evaluate_function(G);
    replace_f_params(Phi,xk_updated,0);  
    TMatrix Phi_k = evaluate_function(Phi);
       
  
    TMatrix P_projected = EKF_P_projection(P_updated,Phi_k,G_k,Q);
    	
    for(nat j = 0; j < rows(out_matrix);j++){ 
      if (j < rows(y)){
        insert_data(out_matrix,j,0,get_data(y_pred,j,0));
      } 
      else{
        insert_data(out_matrix,j,0,get_data(xk_updated,j-rows(y),0));
      }
    }
 /*   
    printf("x_projected: \n");
    display_matrix(x_projected);
*/
    delete_matrix(xk_c);
    delete_matrix(P_c);
    delete_matrix(C_k);
    delete_matrix(U_k);
    delete_matrix(K);
    delete_matrix(yk_m);
    delete_matrix(y_pred);
    delete_matrix(xk_updated);
    delete_matrix(P_updated);
    delete_matrix(G_k);
    delete_matrix(Phi_k);
        
    P_c = P_projected;
    xk_c = x_projected;
    
    replace_f_params(C,xk_c,1);
    replace_f_params(U,xk_c,1);
    replace_f_params(yk,xk_c,1);



    for (nat i = 1; i < cols(y); i++) {
      
     
      C_k = evaluate_function(C);
      U_k = evaluate_function(U);
      /*
      printf("C: \n");
      display_matrix(C_k);
     
      printf("P: \n");
      display_matrix(P_c);
*/
      K = EKF_gain(P_c,C_k,R,U_k);
   /*  
      printf("K: \n");
      display_matrix(K);
     */ 
      y_pred = evaluate_function(yk);
      yk_m = sub_column(y,i);

      xk_updated = EKF_update_state(xk_c,K,yk_m,y_pred); 
      P_updated = update_P(I,P_c,K,C_k);
    

      replace_f_params(xk,xk_updated,i);
      x_projected = evaluate_function(xk);
      
      replace_f_params(G,xk_updated,i);
      G_k = evaluate_function(G);

      replace_f_params(Phi,xk_updated,i);  
      Phi_k = evaluate_function(Phi);

      
      P_projected = EKF_P_projection(P_updated,Phi_k,G_k,Q);
	
      for(nat j = 0; j < rows(out_matrix);j++){ 
        if (j < rows(y)){
          insert_data(out_matrix,j,i,get_data(y_pred,j,0));
        } 
        else{
          insert_data(out_matrix,j,i,get_data(xk_updated,j-rows(y),0));
        }
      }



      delete_matrix(xk_c);
      delete_matrix(P_c);
      delete_matrix(C_k);
      delete_matrix(U_k);
      delete_matrix(K);
      delete_matrix(yk_m);
      delete_matrix(y_pred);
      delete_matrix(xk_updated);
      delete_matrix(P_updated);
      delete_matrix(G_k);
      delete_matrix(Phi_k);
        
      P_c = P_projected;
      xk_c = x_projected;
    
      replace_f_params(C,xk_c,i+1);
      replace_f_params(U,xk_c,i+1);
      replace_f_params(yk,xk_c,i+1);

     
    }
    delete_matrix(P_c);
    delete_matrix(xk_c);
    return out_matrix;

}


