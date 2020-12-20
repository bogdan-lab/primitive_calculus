#include "nonlinear_equations.h"



double func1(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += x[i];
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
    function_t sys[] = {func1, func2, func3};
    double sol[] = {1, 5, 3};
    uint iter_num;
    newton_method(sol, sys, 3, 1e-9, 1e-6, 100, &iter_num);

    printf("Soluiton after %u iterations = ", iter_num);
    for(uint i=0; i<3; i++){
        printf("%.6e ; ", sol[i]);
    }
    printf("\n");
    return 0;
}


