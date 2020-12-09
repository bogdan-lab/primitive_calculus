#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "newton_1D.h"


double bisection(double (*func)(double), const double x_left, const double x_right,
                 const double precision, const uint max_iter, uint* iter_num){
    #ifdef VERBOSE
    printf("**** Bisection method ****\n")
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
        if(idx>=max_iter) break;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx_left = %.6e\tx_right = %.6e\n",
                idx, mid_val, x_left, x_right);
        #endif
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x_mid;
}
