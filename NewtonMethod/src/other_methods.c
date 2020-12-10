﻿#include <stdlib.h>
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
    //Convergence better than linear (~1.5)
    //there are some initial points which cause divergence
    #ifdef VERBOSE
    printf("**** Secant method ****\n");
    #endif
    double val_0 = func(x_0);
    double val_1 = func(x_1);
    double x = (x_0*val_1 - x_1*val_0)/(val_1 - val_0);
    double x_prev = x_1;
    double val_prev = val_1;
    double val  = func(x);
    uint idx = 0;
    while(fabs(val)>precision){
        double correction = (x_prev*val - x*val_prev)/(val - val_prev) - x;
        x_prev = x;
        val_prev = val;
        x = x + correction;
        val = func(x + correction);
        idx++;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e \tcorrection = %.6e\n",
                idx, val, x, correction);
        #endif
        if(idx>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x;
}

void swap(double* lhs, double* rhs){
    double tmp = *lhs;
    *lhs = *rhs;
    *rhs = tmp;
}

double modified_secant_method(double (*func)(double), const double x_0, const double x_1,
                 const double precision, const uint max_iter, uint* iter_num){
    //Trying to implement Newton-like modification
    //SIGNIFICANT INCREASE COMPARE TO CLASSICAL METHOD AND ANOTHER MODIFICATION!!
    //Converges in all cases where classical method does not
    //Also prevents nan
    #ifdef VERBOSE
    printf("**** Modified secant method ****\n");
    #endif
    double val_prev = func(x_0);
    double val = func(x_1);
    double x = x_1;
    double x_prev = x_0;
    if(fabs(val_prev)<fabs(val)){
        swap(&val_prev, &val);
        swap(&x_prev, &x);
    }
    uint idx = 0;
    while(fabs(val)>precision){
        double correction = (x_prev*val - x*val_prev)/(val - val_prev) - x;
        double new_val = func(x+correction);
        while(fabs(new_val)>fabs(val_prev)){
            correction /= 2;
            new_val = func(x+correction);
        }
        if(fabs(new_val)>fabs(val)){
            x_prev = x+correction;
            val_prev = new_val;
        }
        else{
            x_prev = x;
            x = x+correction;
            val_prev = val;
            val = new_val;
        }
        idx++;
        #ifdef VERBOSE
        printf("iteration[%u]\nval_prev = %.6e\tval = %.6e\tx_prev = %.6e \tx = %.6e\n",
                idx, val_prev, val, x_prev, x);
        #endif
        if(idx>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x;
}


double chord_secant_method(double (*func)(double), const double x_0, const double x_1,
                       const double precision, const uint max_iter, uint* iter_num){
    //Seems like in case of chorda method requires well balance of the solution interval
    // if one point is close to the soltion and another - far, than the one which is far away
    //will be moving to towards the solution but angle change in each step will be very small
    // --> bad interval balance - A LOT of chorda iterations!
    //Method converges in some tests where classical was not (in 100 iterations)
    //However there still were some places where it also required more than 100 iterations
    //in some of them classical required less...
    #ifdef VERBOSE
    printf("**** Chord secant method ****\n");
    #endif
    double val_0 = func(x_0);
    double val_1 = func(x_1);
    double x = (x_0*val_1 - x_1*val_0)/(val_1 - val_0);
    double val = func(x);
    double x_prev = x_1;
    double val_prev = val_1;
    uint idx = 0u;
    while(fabs(val)>precision){
        double tmp_x = x;
        x = (x_prev*val - x*val_prev)/(val - val_prev);
        if(val_prev*val<0){
            //solution is between --> do chorda iteration otherwise we will be thrown away
            double new_val = func(x);
            if(val_prev*new_val>0){
                x_prev = tmp_x;
                val_prev = val;
                val = new_val;
            }else {
                val = new_val;
            }
        }else{
            //cannot do chorda iteration --> dobasic secant step
            x_prev = tmp_x;
            val_prev = val;
            val = func(x);
        }
        idx++;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx_prev = %.6e \tx = %.6e\n",
                idx, val, x_prev, x);
        #endif
        if(idx>=max_iter) break;
    }
    if(iter_num != NULL) *iter_num = idx;
    return x;
}
