#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "newton_1D.h"


double simple_func_1(const double x){
    return x+2-exp(x);
    //x1 = -1.841405660; x2 = 1.146193221
}

double simple_func_1_deriv(const double x){
    return 1-exp(x);
}

double simple_func_1_deriv2(const double x){
    return -exp(x);
}

void test_analit(double (*func)(double), double (*deriv)(double),
                 double (*deriv_sec)(double), const uint test_num,
      const double guess_left_bound, const double guess_right_bound){
    for(uint i=0; i<test_num; i++){
        double x = guess_left_bound +
                (guess_right_bound-guess_left_bound)*(1.0*rand()/RAND_MAX);
        double precision = 1e-9;
        uint max_iter = 100;
        printf("TASK: x = %.6e\terror = %.6e\tmax_iter = %u\n",
               x, precision, max_iter);
        uint count = 0;
        double res = third_order_newton_anal(func, x, precision, deriv, deriv_sec,
                                             max_iter, &count);
        printf("CLASSIC: \t root = %.9lf;\tF(root) = %.3e ;\titer_num = %u\n",
               res, func(res), count);
        uint count_mod = 0;
        double res_mod = third_order_newton_modified_anal(func,x, precision,
                   deriv, deriv_sec, max_iter, &count_mod);
        printf("MODIFIED: \t root = %.9lf;\tF(root) = %.3e ;\titer_num = %u\n",
               res, func(res_mod), count_mod);
    }
}


double bad_func_1(const double x){
    if(x==0.0){return 0.0;}
    return x + x*x*sin(2.0/x);
}

void test_calculus(double (*func)(double), const uint test_num,
      const double guess_left_bound, const double guess_right_bound){
    for(uint i=0; i<test_num; i++){
        double x = guess_left_bound +
                (guess_right_bound-guess_left_bound)*(1.0*rand()/RAND_MAX);
        double precision = 1e-9;
        double deriv_delta = 1e-5;
        uint max_iter = 100;
        printf("TASK: x = %.6e\terror = %.6e\tmax_iter = %u\n",
               x, precision, max_iter);
        uint count = 0;
        double res = third_order_newton(func, x, precision, MID_DERIVATIVE,
                     deriv_delta, max_iter, &count);
        printf("CLASSIC: \t root = %.9lf;\tF(root) = %.3e ;\titer_num = %u\n",
               res, func(res), count);
        uint count_mod = 0;
        double res_mod = third_order_newton_modified(func, x, precision, MID_DERIVATIVE,
                     deriv_delta, max_iter, &count_mod);
        printf("MODIFIED: \t root = %.9lf;\tF(root) = %.3e ;\titer_num = %u\n",
               res, func(res_mod), count_mod);
    }
}

double bad_func_2(const double x){
    return x*x;
}



int main(){
    //test_calculus(bad_func_1, 10, -10.0, 10.0);
    //test_calculus(bad_func_2, 10, -10.0, 10.0);
    test_analit(simple_func_1, simple_func_1_deriv, simple_func_1_deriv2,
                10, -10.0, 10.0);
    return 0;
}
