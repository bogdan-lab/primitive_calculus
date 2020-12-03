#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "newton_1D.h"

double simple_func_1(const double x){
    return x+2+exp(x);
    //x1 = -1.841405660; x2 = 1.146193221
}

double simple_func_1_deriv(const double x){
    return 1+exp(x);
}

int main(){
    double x0 = -2.0;
    double precision = 1e-9;
    uint max_iter = 100;
    uint pure_count=0;
    double pure_res = pure_anal_newton_1D(simple_func_1, x0, precision,
                                          simple_func_1_deriv, max_iter,
                                          &pure_count);
    printf("PURE:\t root = %.9lf;\titer_num = %u\n", pure_res, pure_count);
    uint modif_count = 0;
    double modif_res = modified_anal_newton_1D(simple_func_1, x0, precision,
                                          simple_func_1_deriv, max_iter,
                                          &modif_count);
    printf("MODIFIED:\t root = %.9lf;\titer_num = %u\n", modif_res, modif_count);
    return 0;
}
