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
double power(const double x, const uint n){
    double res = 1;
    for(uint i=0; i<n; i++){
        res *= x;
    }
    return res;
}

//FEDORENKO example
/*
double func_f(const double* x, const uint count){
    if(count!=2) exit(1);
    return power(x[0], 5)+ power(x[1], 4) - 2;
}

double func_phi(const double* x, const uint count){
    if(count!=2) exit(1);
    return power(x[0]-2, 3) + power(x[1]-2, 3) + 16;
}
*/

//Ryabenkiy examples
/*
double func_f(const double* x, const uint count){
    if(count!=2) exit(1);
    return sin(x[0]) - x[1] - 1.30;
}

double func_phi(const double* x, const uint count){
    if(count!=2) exit(1);
    return cos(x[1]) - x[0] + 0.84;
}
*/

double func_f(const double* x, const uint count){
    if(count!=2) exit(1);
    return x[0]*x[0] + 4*x[1]*x[1] - 1;
}

double func_phi(const double* x, const uint count){
    if(count!=2) exit(1);
    return power(x[0], 4) + power(x[1], 4) - 0.5;
}



int main(){
    srand(42u);
    const uint eq_num = 2;
    function_t sys[] = {func_f, func_phi};
    //function_t sys[] = {test_linear, test_linear2, test_linear3};
    double sol[eq_num];
    const uint max_iter = 900;
    for(uint i=0; i<1000; i++){
        sol[0] = 500 - (1000.0*rand())/RAND_MAX;
        sol[1] = 500 - (1000.0*rand())/RAND_MAX;
        uint iter_num;
        newton_method(sol, sys, eq_num, 1e-9, 1e-6, max_iter, &iter_num);
        //newton_method_normalized(sol, sys, eq_num, 1e-9, 1e-6, 1000, &iter_num);
        //newton_method_classic(sol, sys, eq_num, 1e-9, 1e-6, 1000, &iter_num);
        if(iter_num>=max_iter){
            printf("ITER NUM = %u\n", iter_num);
            print_result(sol, sys, eq_num);
        }
    }
    return 0;
}


