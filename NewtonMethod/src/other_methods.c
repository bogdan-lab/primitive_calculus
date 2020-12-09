#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "newton_1D.h"


double bisection(double (*func)(double), const double x_left, const double x_right,
                 const double precision, const uint max_iter, uint* iter_num){
    //VERY BAD METHOD it may always converge but it
    //- requires function to change its sign -- huge restriction
    //and with it we need to give boundaries where this sign was changed byt the solution we are looking for
    //Also linear convergence
    #ifdef VERBOSE
    printf("**** Bisection method ****\n");
    #endif
    if(x_left>=x_right){
        fprintf(stderr, "Incorrecct initial interval (%.6e; %.6e)\n", x_left, x_right);
        exit(EXIT_FAILURE);
    }
    double left_val = func(x_left);
    double right_val = func(x_right);
    if(left_val*right_val>0) {
        fprintf(stderr, "Interval (%.6e; %.6e) does not contain solution\n",
                x_left, x_right);
        exit(EXIT_FAILURE);
    }
    if(fabs(left_val)<precision){
        *iter_num = 0;
        return x_left;
    }
    if(fabs(right_val)<precision){
        *iter_num = 0;
        return x_right;
    }
    double length = x_right-x_left;
    double x_start = x_left;
    double x_mid = x_start + 0.5*length;
    double mid_val = func(x_mid);
    uint idx = 0;
    while(fabs(mid_val)>precision){
        length /=2;
        if(mid_val*right_val>0){
            right_val = mid_val;
        } else {
            x_start = x_mid;
        }
        x_mid = x_start + 0.5*length;
        mid_val = func(x_mid);
        idx++;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx_left = %.6e\tx_right = %.6e\n",
                idx, mid_val, x_left, x_right);
        #endif
        if(idx>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x_mid;
}


double secant_method(double (*func)(double), const double x_0, const double x_1,
                 const double precision, const uint max_iter, uint* iter_num){
    #ifdef VERBOSE
    printf("**** Secant method ****\n");
    #endif
    double val_0 = func(x_0);
    double val_1 = func(x_1);
    double x = (x_0*val_1 - x_1*val_0)/(val_1 - val_0);
    double x_prev = fabs(x-x_0)>fabs(x-x_1) ? x_1 : x_0;
    double val_prev = fabs(x-x_0)>fabs(x-x_1) ? val_1 : val_0;
    double val  = func(x);
    uint idx = 0;
    while(fabs(val)>precision){
        double tmp_x = x;
        x = (x_prev*val - x*val_prev)/(val - val_prev);
        x_prev = tmp_x;
        val_prev = val;
        val = func(x);
        idx++;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e \tcorrection = %.6e\n",
                idx, val, x, x-x_prev);
        #endif
        if(idx>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = idx;
    return idx;
}

