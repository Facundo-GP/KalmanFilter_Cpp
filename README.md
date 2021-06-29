
# Kalman filter and Extended Kalman FIter implemented on C++


# Modules description

The objective of this project is to implement a Kalman filter on C++, for this, the work was splited in two modules.


Information about what all the functions implemented do is on .h header files of next modules.


## matrix.cpp

Contains the matrix TAD implementation named TMatrix with all the operations needed to implement the kalman filter.

### kalman.cpp

Contains the Kalman_func TAD, who give a easy way to implement matrix and state update funcionts for EKF. Examples of use are shown on Examples section.

Also it have the implementation of the Kalman filter based on TMatrix types, this module includes funcitons for calculate kalman gain, state projection,P projection, etc for both, Kalman Filter and EKF.

Its also have the implementation of Kalman Filter and EKF itselfs, on the functions Kalman_filter and EKF respectly.


# Examples of use

## Linear model


Lets suposse that we want to estmate the coefficents of a linear model, this can be done by a kalman filter expresing the linear model in a state representation as follows:


EQUATIONS


So all what we need to specify to the function is the values of the diferents matrix, who is shown in the following code:




Finally we call the Kalman Filter function and obtain the output signal and outputs states as shown in the next code




You can vizualize the params and signal predicted by the filter saving the result matrix on a file:



Finally the slope found by the filter is 0.066 for the and here are the evolutions of signal and params over iter steps:


![Demostration](https://drive.google.com/uc?export=view&id=1VkvEppbU7yFRZkLcvo4bKnHUCGSAlVBn)


![Demostration](https://drive.google.com/uc?export=view&id=1i4luAAebZvpbaTyTfpqytUm6aHwt1RYw)

## Sinusoidal wave


![Demostration](https://drive.google.com/uc?export=view&id=11ON2iiE6WQNzhqprtJ-JkWjEdEfUYc5r)


![Demostration](https://drive.google.com/uc?export=view&id=1B-tGxxEEsHJ9ylYBNJVu3I51Z9axKbsS)

## Second order system



![Demostration](https://drive.google.com/uc?export=view&id=100QZJVQxt9qAOwx-DNzQRPVl6DTVe1I7)


![Demostration](https://drive.google.com/uc?export=view&id=1QA1oxbRlMf7DLzVpkw8T2DyDHiHZSm3u)
