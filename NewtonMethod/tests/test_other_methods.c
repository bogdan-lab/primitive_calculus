#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "newton_1D.h"

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

void bisection_test(double (*func)(double), const uint test_num,
      const double left_guess, const double right_guess, const double guess_affinity){
    for(uint i=0; i<test_num; i++){
        double x_left = left_guess+ (2.0*rand()/RAND_MAX -1.0)*guess_affinity;
        double x_right = right_guess+ (2.0*rand()/RAND_MAX -1.0)*guess_affinity;
        double precision = 1e-9;
        uint max_iter = 100;
        printf("TASK: xl = %.6e\tx_r = %.6e\terror = %.6e\tmax_iter = %u\n",
               x_left, x_right, precision, max_iter);
        uint count = 0;
        double res = bisection(func, x_left, x_right, precision, max_iter, &count);
        printf("RES:\t root = %.9lf;\tF(root) = %.3e ;\titer_num = %u\n",
               res, func(res), count);
    }
}


int main(){
    const uint test_num = 100;
    double left_guess = 0.0;
    double righ_guess = M_PI;
    double guess_affinity = 0.5*M_PI;
    //bisection_test(simple_func_1, test_num, left_guess, righ_guess, guess_affinity);
    left_guess = -3.0;
    righ_guess = -0.5;
    guess_affinity = 1.0;
    //bisection_test(simple_func_2, test_num, left_guess, righ_guess, guess_affinity);
    left_guess = -5;
    righ_guess = 5;
    guess_affinity = 2;
    bisection_test(bad_func_1, test_num, left_guess, righ_guess, guess_affinity);
    return 0;
}
