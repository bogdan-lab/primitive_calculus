#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nonlinear_equation_1D.h"

double get_derivative(double (*func)(), const double x0, const double delta,
                      const DerivType der_type){
    switch (der_type) {
    case LEFT_DERIVATIVE:
        return (func(x0) - func(x0-delta))/delta;
    case RIGHT_DERIVATIVE:
        return (func(x0 + delta)-func(x0))/delta;
    case MID_DERIVATIVE:;
        return (func(x0 + delta) - func(x0 - delta))/(2*delta);
    default:
        fprintf(stderr, "Incorrect derivative type flag\n");
        exit(EXIT_FAILURE);
    }
}

double pure_calc_newton_1D(double (*func)(double), const double x0, const double precision,
   const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_made){
    #ifdef VERBOSE
    printf("**** Pure Numerical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(fabs(defect)>precision){
        double derivative = get_derivative(func, x, der_delta, der_type);
        double correction = defect/derivative;
        x -= correction;
        idx++;
        if(idx>=max_iter){
            break;
        }
        defect = func(x);
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, defect, x, correction);
        #endif
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;
}

double pure_anal_newton_1D(double (*func)(double), const double x0, const double precision,
                       double (*deriv)(double), const uint max_iter, uint* iter_made){
    #ifdef VERBOSE
    printf("**** Pure Analitical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(fabs(defect)>precision){
        double correction = defect/deriv(x);
        x -= correction;
        idx++;
        if(idx>=max_iter){
            break;
        }
        defect = func(x);
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, defect, x, correction);
        #endif
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;
}



//modified newton with dividng by two
double modified_calc_newton_1D(double (*func)(double), const double x0, const double precision,
    const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_made){
    #ifdef VERBOSE
    printf("**** Modified Numerical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(fabs(defect)>precision){
        double derivative = get_derivative(func, x, der_delta, der_type);
        double correction = -defect/derivative;
        double new_val = func(x + correction);
        while(fabs(new_val)>fabs(defect)){
            correction /= 2;
            new_val = func(x+correction);
        }
        x += correction;
        idx++;
        if(idx>=max_iter){
            break;
        }
        defect = new_val;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, defect, x, correction);
        #endif
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;

}

double modified_anal_newton_1D(double (*func)(double), const double x0, const double precision,
                       double (*deriv)(double), const uint max_iter, uint* iter_made){
    #ifdef VERBOSE
    printf("**** Modified Analitical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(fabs(defect)>precision){
        double correction = -defect/deriv(x);
        double new_val = func(x+correction);
        while(fabs(new_val)>fabs(defect)){
            correction /= 2;
            new_val = func(x+correction);
        }
        x += correction;
        idx++;
        if(idx>=max_iter){
            break;
        }
        defect = new_val;
        #ifdef VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\tcorrection = %.6e\n",
               idx, defect, x, correction);
        #endif
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;

}

