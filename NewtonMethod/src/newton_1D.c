#include <stdlib.h>
#include <stdio.h>
#include "newton_1D.h"

double pure_calc_newton_1D(double (*func)(double), const double x0, const double precision,
   const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_made){
    #ifdef NEWTON_VERBOSE
    printf("**** Pure Numerical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(abs(defect)>precision){
        double derivative = 0;
        switch (der_type) {
        case LEFT_DERIVATIVE:
            derivative = (defect - func(x-der_delta))/der_delta;
            break;
        case RIGHT_DERIVATIVE:
            derivative = (func(x+der_delta)-defect)/der_delta;
            break;
        case MID_DERIVATIVE:;
            derivative = (func(x+der_delta) - func(x-der_delta))/der_delta;
            break;
        default:
            fprintf(stderr, "Incorrect derivative type flag\n");
            exit(EXIT_FAILURE);
        }
        x -= defect/derivative;
        idx++;
        #ifdef NEWTON_VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\n", idx, defect);
        #endif
        if(idx>=max_iter){
            break;
        }
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;
}

double pure_anal_newton_1D(double (*func)(double), const double x0, const double precision,
                       double (*deriv)(double), const uint max_iter, uint* iter_made){
    #ifdef NEWTON_VERBOSE
    printf("**** Pure Analitical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(abs(defect)>precision){
        x -= defect/deriv(x);
        idx++;
        #ifdef NEWTON_VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\n", idx, defect);
        #endif
        if(idx>=max_iter){
            break;
        }
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;
}



//modified newton with dividng by two
double modified_calc_newton_1D(double (*func)(double), const double x0, const double precision,
    const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_made){
    #ifdef NEWTON_VERBOSE
    printf("**** Modified Numerical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(abs(defect)>precision){
        double derivative = 0;
        switch (der_type) {
        case LEFT_DERIVATIVE:
            derivative = (defect - func(x-der_delta))/der_delta;
            break;
        case RIGHT_DERIVATIVE:
            derivative = (func(x+der_delta)-defect)/der_delta;
            break;
        case MID_DERIVATIVE:;
            derivative = (func(x+der_delta) - func(x-der_delta))/der_delta;
            break;
        default:
            fprintf(stderr, "Incorrect derivative type flag\n");
            exit(EXIT_FAILURE);
        }
        double correction = defect/derivative;
        double new_val = x - correction;
        while(func(new_val)>func(x)){
            correction /= 2;
            new_val += correction;
        }
        x = new_val;
        idx++;
        #ifdef NEWTON_VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\n", idx, defect);
        #endif
        if(idx>=max_iter){
            break;
        }
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;

}

double modified_anal_newton_1D(double (*func)(double), const double x0, const double precision,
                       double (*deriv)(double), const uint max_iter, uint* iter_made){
    #ifdef NEWTON_VERBOSE
    printf("**** Modified Analitical Newton 1D ****\n");
    #endif
    double x = x0;
    double defect = func(x);
    uint idx = 0;
    while(abs(defect)>precision){
        double correction = defect/deriv(x);
        double new_val = x - correction;
        while(func(new_val)>func(x)){
            correction /= 2;
            new_val += correction;
        }
        x = new_val;
        idx++;
        #ifdef NEWTON_VERBOSE
        printf("iteration[%u]\ndefect = %.6e\tx = %.6e\n", idx, defect);
        #endif
        if(idx>=max_iter){
            break;
        }
    }
    if(iter_made!=NULL) *iter_made = idx;
    return x;

}

//other methods: sekushchie, bisection..
// find all roots by newton deleting other

//normalized newton for systems
