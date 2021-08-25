



# Kalman Filter and Extended Kalman Filter implemented on C++.

 <br>

# Modules description
 

The objective of this project is to implement a Kalman filter on C++, for this, the work was split in two modules.

Information about what all the functions implemented do is on .h header files of next modules.

For compilation you can build a make file or if you using g++ compile with g++  source/matrix.cpp source/kalman.cpp your_program.cpp -o executable


 <br>

## examples directory

Contains the .cpp files of all examples named on "Examples of use" section.



 <br>

## source directory

Contains .cpp files.

 <br>

### matrix.cpp
Contains the matrix TAD implementation named TMatrix with all the operations needed to implement the kalman filter.

 <br>

### kalman.cpp

Contains the Kalman_func TAD, who give an easy way to implement matrix and state update functions for EKF. Examples of use are shown on Examples section.

Also, it has the implementation of the Kalman filter based on TMatrix types, this module includes functions to calculate kalman gain, state projection, P projection, and others for both, Kalman Filter and EKF.

It's also had the implementation of Kalman Filter and EKF itself, on the functions Kalman_filter and EKF respectively.

 <br>

## include directory

Contains header files.
 
 <br>
  <br>

# Examples of use

For the state space representations of variables that Kalman filter use, the follow notation will be used:

 <br>


![E1](https://drive.google.com/uc?export=view&id=1gfH2q1sz3MSRzHJMVrn0MFkCXUuhvYs9)


 <br>

Autocorrelation matrices of vk and wk will be denoted R and Q respectly.

Finally, for EKF model, we will use:

 <br>

![E2](https://drive.google.com/uc?export=view&id=15pWBq5VsvSmy7Pgf-uxLUMu0MPhvBsh3)

 <br>

And the jacobian matrices:

 <br>

 
 ![E3](https://drive.google.com/uc?export=view&id=1wCIzRkDw5uPmhP38kB2-Fv1Vkj2ib6Xd)
 
  <br>
  
  
 ![E4](https://drive.google.com/uc?export=view&id=19G1NzqCd_qH48Rp42Qi7nUs0IwNLb3HV)
 
 <br>

## Linear model

 <br>

Let's suppose that you want to estimate the coefficients of a linear model, this can be done by a kalman filter expressing the linear model in a state representation as follows:

 <br>

![E5](https://drive.google.com/uc?export=view&id=1z9Usgy4UmrqIZ6pFK88gj23xmxpKjduJ)
 


 <br>

In this case,[r,p] transpose is your state vector and vk white Gaussian noise with  mean 0 and variance 9.
This model represents a linear equation. The objective is then estimating pk state who is supposed constant, and represents the slope of the line.

So all you need to specify to the function is the values of C, Phi ,R and Q matrices and the prior of the state and P error matrix, who is shown in the following code:


 <br>

![Code_linear1](https://drive.google.com/uc?export=view&id=194Y5wn_monaPXIfNAiVB6VNTAcgli9m-)


 <br>

Where nat type is an unsigned int.
 
Finally,call the Kalman Filter function and obtain the output signal and outputs states as shown in the next code

<br>
<br>

![Code_linear2](https://drive.google.com/uc?export=view&id=1i6ioGYe5l978hsp2w1SwodKS-bepQksR)

 <br>
  <br>

Where "I" matrix is the identity of size (num_states , num_states) and "y" the TMatrix with the observed signal.

The resulting matrices of these functions are TMatrix of dimension (n_outputs, num_iterations) for signal and (num_states, num_iterations) for states.

You can visualize the params and signal predicted by the filter, saving the result matrix on a file with get_data function from matrix.h file and plotting it on a software you want.

Finally, the slope found by the filter is 0.066 for the and here are the evolutions of signal and params over iteration steps:
:

<br>
<br>

![Signal_linear](https://drive.google.com/uc?export=view&id=1VkvEppbU7yFRZkLcvo4bKnHUCGSAlVBn)

 <br>
 <br> 
 <br>

![States_linear](https://drive.google.com/uc?export=view&id=1i4luAAebZvpbaTyTfpqytUm6aHwt1RYw)


<br>
<br> 
<br>


## Sinusoidal wave

Suppose now you want to estimate the amplitude and phase of a sinusoidal wave that you know the frequency who measures are  with white Gaussian noise with mean 0 and variance 1 added.

In this case it's necessary to use the EKF and specify the C ,Phi, G, U, state and output evolution equation.

This problem can be represented in state space as follows:

<br>
<br>


![E6](https://drive.google.com/uc?export=view&id=1Uk-L4xLlANJdpdwTiXkbZblj1fRVruLs)


<br>
<br>

So, C matrix is:

![E7](https://drive.google.com/uc?export=view&id=1LWJKm6TfgfJNtnzSgpGe6cE9lThnyQnL)

Phi and U are the identity matrix, Q is the zero matrix and observed output is:

![E8](https://drive.google.com/uc?export=view&id=1kJlzZeZTx3t82lzTDdqoY8PdGERw6qpU)

Since you are working with EKF, matrices and evolution functions may be of TEKF_func type. Assuming a frequency of pi/5  and a sample_rate of 0.25 you can write it as follows:


<br>
<br>

![Demostration](https://drive.google.com/uc?export=view&id=19j68dAEUf7VzTEpQdowLjsJfhhB68P4A)


Where x_prior is the initial state, 0 is the start iteration and others arguments are the evolution function of the different matrices, here is an example of declaration of C_k and y_k for this problem:



![Demostration](https://drive.google.com/uc?export=view&id=1vl7CtJ6qL5zd3cW3Rwui9l7gh6N3Lszt)


<br>
<br>

Where state_params the state vector and iter_param the iteration, inside EKF function these parameters are updated on each call to the functions with the predicted values.

Finally, you can call EKF function as follows:

<br>
<br>


![Demostration](https://drive.google.com/uc?export=view&id=1q_eLREeKRmCsFEUBmvgEUV7cBOhQs922)


<br>
<br>


Where I, P_prior, Q, R, y and p_prior are TMatrix type, since them don't change over iterations.

As on linear model, I is the identity of size (n_states, n_states).

It's important to make a difference, while "y" parameter is a TMatrix variable containing the measured function, "yk" is a TEKF_func variable who haves the observed output evolution function y_k inside.


You can save the result like on linear model example and plot it, in this case EKF estimates 1.907 for amplitude and 1.09 for phase, while true values were 2 for amplitude and pi/3 for phase.


<br>
<br> 
<br>

![Demostration](https://drive.google.com/uc?export=view&id=11ON2iiE6WQNzhqprtJ-JkWjEdEfUYc5r)



<br>
<br> 
<br>

![Demostration](https://drive.google.com/uc?export=view&id=1B-tGxxEEsHJ9ylYBNJVu3I51Z9axKbsS)



<br>
<br> 
<br>


## Second order system

For  last example, suppose that you want to estimate the ξ , A and Wn of this second order system assuming that ξ << 1 and 0 < ξ < 1 :

<br>
<br>


![E9](https://drive.google.com/uc?export=view&id=1qQU3moArH6QHDzSoxptb8TtrJZx-4e_A)


<br>
<br>


For this, you measure the step response which assuming ξ << 1 can be expressed in a state space model like this:



<br>
<br>


![E10](https://drive.google.com/uc?export=view&id=1GpQTFV7AYgiWgPLAVWH_P6xv9VZ-n3_x)



<br>
<br>

Where vk is white Gaussian noise of mean 0 and variance 4.

In this model, the state vector is [Ak,ξk,Wnk] transpose. Since you know that $\xi$ is between 0 and 1, you can put it as a restriction on the state evolution function as follows:




<br>
<br>


![Demostration](https://drive.google.com/uc?export=view&id=1Woe5LU7wszfX3nptYXpZm2I4cMyokZGa)

<br>
<br>



The rest of the implementation is equal to the sinusoidal example. Values used to generate the signal were 6 for A, 0.1 for ξ and 15 for Wn, while the estimated by EKF were 8.04 for A, 0.107 for ξ and 14.83 for Wn.

<br>
<br>


![Demostration](https://drive.google.com/uc?export=view&id=100QZJVQxt9qAOwx-DNzQRPVl6DTVe1I7)



<br>
<br> 
<br>

![Demostration](https://drive.google.com/uc?export=view&id=1QA1oxbRlMf7DLzVpkw8T2DyDHiHZSm3u)



