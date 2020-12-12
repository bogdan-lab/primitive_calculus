#include <stdio.h>

#include "newton_1D.h"


double get_second_derivative(double (*func)(double), const double x,
                             const double der_delta){
    return (func(x+der_delta) - 2*func(x) + func(x-der_delta))/(der_delta*der_delta);
}

double third_order_newton(double (*func)(double), const double x0, const double precision,
  const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_num){
    #ifdef VERBOSE
    printf("**** Third order Newton method Calculus****\n")
    #endif
    double x = x0;
    double val = func(x0);
    uint idx = 0;
    while(fabs(val)>precision){
        double der = get_derivative(func, x, der_delta, der_type);
        double der_sec = get_second_derivative(func, x, der_delta);
        double correction = -val/der - 0.5*der_sec*val*val/(der*der*der);
        x+=correction;
        idx++;
        if(idx>=max_iter) break;
        val = func(x);
        #ifdef VERBOSE
        printf("iteration[%u]\tdefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, val, x, correction);
        #endif
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x;
}

double third_order_newton_anal(double (*func)(double), const double x0,
  const double precision, double (*deriv)(double), double (*deriv_sec)(double),
       const uint max_iter, uint* iter_num){
    #ifdef VERBOSE
    printf("**** Third order Newton method Analitical ****\n")
    #endif
    double x = x0;
    double val = func(x0);
    uint idx = 0;
    while(fabs(val)>precision){
        double der = deriv(x);
        double der_sec = deriv_sec(x);
        double correction = -val/der - 0.5*der_sec*val*val/(der*der*der);
        x+=correction;
        idx++;
        if(idx>=max_iter) break;
        val = func(x);
        #ifdef VERBOSE
        printf("iteration[%u]\tdefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, val, x, correction);
        #endif
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x;
}


double third_order_newton_modified(double (*func)(double), const double x0,
        const double precision, const DerivType der_type, const double der_delta,
           const uint max_iter, uint* iter_num){
    #ifdef VERBOSE
    printf("**** Third order Newton method Calculus modified****\n")
    #endif
    double x = x0;
    double val = func(x0);
    uint idx = 0;
    while(fabs(val)>precision){
        double der = get_derivative(func, x, der_delta, der_type);
        double der_sec = get_second_derivative(func, x, der_delta);
        double correction = -val/der - 0.5*der_sec*val*val/(der*der*der);
        double new_val = func(x+correction);
        while (fabs(new_val)>fabs(val)) {
            correction /= 2;
            new_val = func(x+correction);
        }
        x+=correction;
        idx++;
        if(idx>=max_iter) break;
        val = new_val;
        #ifdef VERBOSE
        printf("iteration[%u]\tdefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, val, x, correction);
        #endif
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x;
}


double third_order_newton_modified_anal(double (*func)(double), const double x0,
     const double precision, double (*deriv)(double), double (*deriv_sec)(double),
           const uint max_iter, uint* iter_num){
    #ifdef VERBOSE
    printf("**** Third order Newton method Analytical modified****\n")
    #endif
    double x = x0;
    double val = func(x0);
    uint idx = 0;
    while(fabs(val)>precision){
        double der = deriv(x);
        double der_sec = deriv_sec(x);
        double correction = -val/der - 0.5*der_sec*val*val/(der*der*der);
        double new_val = func(x+correction);
        while (fabs(new_val)>fabs(val)) {
            correction /= 2;
            new_val = func(x+correction);
        }
        x+=correction;
        idx++;
        if(idx>=max_iter) break;
        val = new_val;
        #ifdef VERBOSE
        printf("iteration[%u]\tdefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, val, x, correction);
        #endif
    }
    if(iter_num!=NULL) *iter_num = idx;
    return x;
}



