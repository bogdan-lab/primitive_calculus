#include <stdio.h>
#include <math.h>

#include "newton_1D.h"


double get_second_derivative(double (*func)(double), const double x,
                             const double der_delta){
    return (func(x+der_delta) - 2*func(x) + func(x-der_delta))/(der_delta*der_delta);
}

double third_order_newton(double (*func)(double), const double x0, const double precision,
  const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_num){
    //SHOULD BE VERY CAREFUL with derivatieves!
    //in classical newton f(x)/f'(x) requires only that f'(x) is decreases more slowely than f(x)
    //Here we have f''(x) in the numerator --> if close to the solution it become infinite we hav problem
    //ALSO IN CASE OF CALCULATIVE DERIVATIVE much more sensitive to it than classical second order newton
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
    //GOOD IF ALL DERIVATIVES NEAR THE SOLUTION ARE OK!
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
    //IN MY POOR TESTS HAD THE SAME RESULTS OR WORSE THAN THOSE OBTAINED BY NONMODIFIED
    //Results become worse in case with crazy infinit changeable second derivative near the solution
    //Funny in that case modification simply clamps iterations close to the point where f'' start diverge
    //And from that point correction on each iteration becomes very small
    //Nonmodified method sometimes converges in such "crazy f''" case because of the "fly effect"
    //it simply is blown away from the bad position by one iteration and starts over and some
    //of those start-overs finishes with success. However such succes is ocasional and should not be relied on!
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



