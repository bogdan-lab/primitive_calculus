﻿#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "nonlinear_equation_1D.h"

double simple_func_1(const double x){
    //Cannot be solved by bisection because it is never less than zero!
    return (sin(x) - 1)*(sin(x)-1);
    //solutions = n*M_PI/2
}

double simple_func_2(const double x){
    return x+2-exp(x);
    //x1 = -1.841405660; x2 = 1.146193221
}

double bad_func_1(const double x){
    if(x==0.0) return 0.0;
    return x + x*x*sin(2.0/x);
}


void secant_test(double (*func)(double), const uint test_num,
      const double left_guess, const double right_guess, const double guess_affinity){
    uint count_bad_basic = 0;
    uint count_bad_modified = 0;
    for(uint i=0; i<test_num; i++){
        double x_0 = left_guess+ (2.0*rand()/RAND_MAX -1.0)*guess_affinity;
        double x_1 = right_guess+ (2.0*rand()/RAND_MAX -1.0)*guess_affinity;
        double precision = 1e-9;
        uint max_iter = 100;
        printf("TASK: x_0 = %.6e\tx_1 = %.6e\terror = %.6e\tmax_iter = %u\n",
               x_0, x_1, precision, max_iter);
        uint count = 0;
        double res = secant_method(func, x_0, x_1, precision, max_iter, &count);
        printf("BASIC:\t root = %.9lf;\tF(root) = %.6e ;\titer_num = %u\n",
               res, func(res), count);
        if(count==max_iter) count_bad_basic++;
        uint mod_count = 0;
        double mod_res = modified_secant_method(func, x_0, x_1, precision, max_iter, &mod_count);
        printf("MODIFIED:\t root = %.9lf;\tF(root) = %.6e ;\titer_num = %u\n",
               mod_res, func(mod_res), mod_count);
        if(mod_count==max_iter) count_bad_modified++;
    }
    printf("Number of bad tests:\nBasic = %u;\n Modified = %u\n",
           count_bad_basic, count_bad_modified);
}


int main(){
    //SECANT METHOD
    //secant_test(simple_func_1, 1000, 0.0, 1.0, 0.5*M_PI);
    secant_test(simple_func_2, 1000, 0.0, 10.0, 1);
    //secant_test(bad_func_1, 1000, -10.0, 10.0, 5.0);
    return 0;
}
