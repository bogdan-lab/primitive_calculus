#include <stdio.h>

#include "stdlib.h"
#include "lapacke.h"
#include "gsl/gsl_cblas.h"

#include "nonlinear_equations.h"

double get_derivative_val(function_t func, const uint pos, const double* point,
                          const uint count, const double der_delta){
    double* lhs_point = calloc(count, sizeof(*lhs_point));
    double* rhs_point = calloc(count, sizeof(*lhs_point));
    memcpy(lhs_point, point, count*sizeof(double));
    memcpy(rhs_point, point, count*sizeof(double));
    lhs_point[pos] += der_delta;
    rhs_point[pos] -= der_delta;
    return (func(lhs_point, count) - func(rhs_point, count))/(2*der_delta);
}

void fill_derivative_matrix(function_t* system, const double* point, const uint count,
                            const double der_delta, double* der_matrix){
    const uint m_size = count*count;
    for(uint i=0; i<m_size; i++){
        uint func_idx = i/count;
        uint arg_idx = i % count;
        der_matrix[i] = get_derivative_val(system[func_idx], arg_idx, point,
                                           count, der_delta);
    }
}

int check_solution(function_t* system, const double* solution, const uint count,
                   const double precision){
    double sum = 0;
    for(uint i=0; i<count ; i++){
        sum += fabs(system[i](solution, count));
    }
    return sum<precision;
}

void fill_func_val(function_t* system, const double* point, const uint count,
                   double* val_vec){
    for(uint i=0; i<count; i++){
        val_vec[i] = system[i](point, count);
    }
}

void inverse_square_matrix(double* matrix, const uint n){
    int ipiv[n];
    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matrix, n, ipiv);
    error_hendler(info);
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, matrix, n, ipiv);
    error_hendler(info);
}

void newton_method(double* solution, function_t* system, const uint count,
                    const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num){
    uint iter = 0;
    while(!check_solution(system, solution, count, precision)){
        double der_matrix[count*count];
        fill_derivative_matrix(system, solution, count, der_delta, der_matrix);
        inverse_square_matrix(der_matrix, count);
        double rhs[count];
        fill_func_val(system, solution, count, rhs);
        double delta_x[count];
        cblas_dgemv(CblasRowMajor, CblasNoTrans, count, count, -1.0, der_matrix,
                    count, rhs, 1, 0.0, delta_x, 1);
        cblas_daxpy(count, 1.0, delta_x, 1, solution, 1);
        iter++;
        if(iter>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = iter;
}

void error_hendler(const int info){
    if(info!=0){
        fprintf(stderr, "Something went wrong!");
        exit(EXIT_FAILURE);
    }
}

