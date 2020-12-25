#include "nonlinear_equations.h"


/*
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

double test_linear(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += x[i];
    }
    return sum;
}

double test_linear2(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += 2*x[i];
    }
    return sum;
}

double test_linear3(const double* x, const uint count){
    double sum = 0;
    for(uint i=0; i<count; i++){
        sum += 3*x[i];
    }
    return sum;
}
*/
//Book example
double power(const double x, const uint n){
    double res = 1;
    for(uint i=0; i<n; i++){
        res *= x;
    }
    return res;
}

double func_f(const double* x, const uint count){
    if(count!=2) exit(1);
    return power(x[0], 5)+ power(x[1], 4) - 2;
}

double func_phi(const double* x, const uint count){
    if(count!=2) exit(1);
    return power(x[0]-2, 3) + power(x[1]-2, 3) + 16;
}

int main(){
    srand(42u);
    const uint eq_num = 2;
    function_t sys[] = {func_f, func_phi};
    //function_t sys[] = {test_linear, test_linear2, test_linear3};
    double sol[eq_num];
    for(uint i=0; i<5; i++){
        sol[0] = 2; //(10.0*rand())/RAND_MAX;
        sol[1] = 3; //(10.0*rand())/RAND_MAX;
        uint iter_num;
        //newton_method(sol, sys, eq_num, 1e-9, 1e-6, 1000, &iter_num);
        newton_method_normalized(sol, sys, eq_num, 1e-9, 1e-6, 1000, &iter_num);
        //newton_method_classic(sol, sys, eq_num, 1e-9, 1e-6, 1000, &iter_num);
        printf("ITER NUM = %u\n", iter_num);
        print_result(sol, sys, eq_num);
    }
    return 0;
}


