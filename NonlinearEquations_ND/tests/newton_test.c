#include "nonlinear_equations.h"



double func1(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += cbrt(x[i]);
    }
    return sum;
}

double func2(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += x[i]*x[i];
    }
    return sum;
}

double func3(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += x[i]*x[i]*x[i];
    }
    return sum;
}


int main(){
    srand(43u);
    const uint eq_num = 3;
    function_t sys[] = {func1, func2, func3};
    double sol[] = {0.1, 0.1, 0.1};
    uint iter_num;
    newton_method(sol, sys, eq_num, 1e-9, 1e-6, 100, &iter_num);
    //newton_method_classic(sol, sys, eq_num, 1e-9, 1e-6, 100, &iter_num);
    printf("ITER NUM = %u\n", iter_num);
    print_result(sol, sys, eq_num);
    return 0;
}


