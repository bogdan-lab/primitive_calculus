#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "newton_1D.h"

double simple_func_1(const double x){
    return x+2-exp(x);
    //x1 = -1.841405660; x2 = 1.146193221
}

double simple_func_1_deriv(const double x){
    return 1-exp(x);
}

void test_simple_func_1(const uint test_num, const double left_guess_bound,
                        const double right_guess_bound){
    for(uint i=0; i<test_num; i++){
        double x0 = left_guess_bound + (right_guess_bound - left_guess_bound)*1.0*rand()/RAND_MAX;
        //double x0 = 7.419549e-4; 	<--Pure gives nan since we are close to deriv = 0
        double precision = 1e-9;
        uint max_iter = 100;
        printf("TASK: x0 = %.6e\terror = %.6e\tmax_iter = %u\n", x0, precision, max_iter);
        uint pure_count=0;
        double pure_res = pure_anal_newton_1D(simple_func_1, x0, precision,
                                              simple_func_1_deriv, max_iter,
                                              &pure_count);
        printf("PURE:\t root = %.9lf;\tF(root) = %.3e;\titer_num = %u\n", pure_res, simple_func_1(pure_res), pure_count);
        uint modif_count = 0;
        double modif_res = modified_anal_newton_1D(simple_func_1, x0, precision,
                                                   simple_func_1_deriv, max_iter,
                                                   &modif_count);
        printf("MODIFIED:\t root = %.9lf\tF(root) = %.3e;\t\titer_num = %u\n", modif_res, simple_func_1(modif_res), modif_count);
        printf("\n-----------------------------------------------------------------------\n");
    }
}

double bad_func_1(const double x){
    //Crazy derivative randomly changes sign of correction with approaching
    //solution --> very problematic for Newton, may be easy for non-derivative method
    //still modifed method is much better than pure
    if(x==0.0){return 0.0;}
    return x + x*x*sin(2.0/x);
}

void test_bad_func_1(const uint test_num, const double left_guess_bound,
                        const double right_guess_bound){
    for(uint i=0; i<test_num; i++){
        double x0 = left_guess_bound + (right_guess_bound - left_guess_bound)*1.0*rand()/RAND_MAX;
        double precision = 1e-9;
        uint max_iter = 1000;
        DerivType deriv_type = MID_DERIVATIVE;
        double deriv_delta = precision;
        printf("TASK: x0 = %.6e\terror = %.6e\tmax_iter = %u\n", x0, precision, max_iter);
        uint pure_count=0;
        double pure_res = pure_calc_newton_1D(bad_func_1, x0, precision,
                                              deriv_type, deriv_delta,
                                              max_iter, &pure_count);
        printf("PURE:\t root = %.9e;\tF(root) = %.3e;\titer_num = %u\n", pure_res, bad_func_1(pure_res), pure_count);
        uint modif_count = 0;
        double modif_res = modified_calc_newton_1D(bad_func_1, x0, precision,
                                                   deriv_type, deriv_delta,
                                                   max_iter, &modif_count);
        printf("MODIFIED:\t root = %.9e\tF(root) = %.3e;\t\titer_num = %u\n", modif_res, bad_func_1(modif_res), modif_count);
        printf("\n-----------------------------------------------------------------------\n");
    }
}


double bad_func_2(const double x){
    //No second order derivative in the root --> slow convergence
    //Do not really feel it
    return x + x*cbrt(x);
}

double bad_func_2_deriv(const double x){
    return 1 + 4.0/3.0*cbrt(x);
}

void test_bad_func_2(const uint test_num, const double left_guess_bound,
                        const double right_guess_bound){
    for(uint i=0; i<test_num; i++){
        double x0 = left_guess_bound + (right_guess_bound - left_guess_bound)*1.0*rand()/RAND_MAX;
        double precision = 1e-9;
        uint max_iter = 100;
        printf("TASK: x0 = %.6e\terror = %.6e\tmax_iter = %u\n", x0, precision, max_iter);
        uint pure_count=0;
        double pure_res = pure_anal_newton_1D(bad_func_2, x0, precision,
                                              bad_func_2_deriv, max_iter,
                                              &pure_count);
        printf("PURE:\t root = %.9lf;\tF(root) = %.3e;\titer_num = %u\n", pure_res, bad_func_2(pure_res), pure_count);
        uint modif_count = 0;
        double modif_res = modified_anal_newton_1D(bad_func_2, x0, precision,
                                                   bad_func_2_deriv, max_iter,
                                                   &modif_count);
        printf("MODIFIED:\t root = %.9lf\tF(root) = %.3e;\t\titer_num = %u\n", modif_res, bad_func_2(modif_res), modif_count);
        printf("\n-----------------------------------------------------------------------\n");
    }
}


double bad_func_3(const double x){
    //If derivative in the root is 0 --> convergence is not squared --> really see increase of iteration number
    //obviously can have noticeable mistake in solution <-- A lot points will give close to zero results
    return x*x;
}

double bad_func_3_deriv(const double x){
    return 2*x;
}

void test_bad_func_3(const uint test_num, const double left_guess_bound,
                        const double right_guess_bound){
    for(uint i=0; i<test_num; i++){
        double x0 = left_guess_bound + (right_guess_bound - left_guess_bound)*1.0*rand()/RAND_MAX;
        double precision = 1e-9;
        uint max_iter = 100;
        printf("TASK: x0 = %.6e\terror = %.6e\tmax_iter = %u\n", x0, precision, max_iter);
        uint pure_count=0;
        double pure_res = pure_anal_newton_1D(bad_func_3, x0, precision,
                                              bad_func_3_deriv, max_iter,
                                              &pure_count);
        printf("PURE:\t root = %.9lf;\tF(root) = %.3e;\titer_num = %u\n", pure_res, bad_func_3(pure_res), pure_count);
        uint modif_count = 0;
        double modif_res = modified_anal_newton_1D(bad_func_3, x0, precision,
                                                   bad_func_3_deriv, max_iter,
                                                   &modif_count);
        printf("MODIFIED:\t root = %.9lf\tF(root) = %.3e;\t\titer_num = %u\n", modif_res, bad_func_3(modif_res), modif_count);
        printf("\n-----------------------------------------------------------------------\n");
    }
}

int main(){
    uint test_num = 10000;
    double left_guess_bound = -10.0;
    double right_guess_bound = 10.0;
    srand((uint)time(NULL));
    //test_simple_func_1(test_num, left_guess_bound, right_guess_bound);
    //test_bad_func_1(test_num, left_guess_bound, right_guess_bound);
    //test_bad_func_2(test_num, left_guess_bound, right_guess_bound);
    test_bad_func_3(test_num, left_guess_bound, right_guess_bound);
    return 0;
}
