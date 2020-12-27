#include <stdio.h>

#include "stdlib.h"
#include "lapacke.h"
#include "gsl/gsl_cblas.h"

#include "nonlinear_equations.h"

double get_derivative_val(const function_t func, const uint pos, const double* point,
                          const uint count, const double der_delta){
    double lhs_point[count];
    double rhs_point[count];
    memcpy(lhs_point, point, count*sizeof(double));
    memcpy(rhs_point, point, count*sizeof(double));
    lhs_point[pos] += der_delta;
    rhs_point[pos] -= der_delta;
    return (func(lhs_point, count) - func(rhs_point, count))/(2*der_delta);
}

void fill_derivative_matrix(const function_t* system, const double* point, const uint count,
                            const double der_delta, double* der_matrix){
    const uint m_size = count*count;
    for(uint i=0; i<m_size; i++){
        uint func_idx = i/count;
        uint arg_idx = i % count;
        der_matrix[i] = get_derivative_val(system[func_idx], arg_idx, point,
                                           count, der_delta);
    }
}

int check_solution(const function_t* system, const double* solution, const uint count,
                   const double precision){
    double sum = 0;
    for(uint i=0; i<count ; i++){
        sum += fabs(system[i](solution, count));
    }
    return sum<precision;
}

void fill_func_val(const function_t* system, const double* point, const uint count,
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

void print_result(const double* solution, const function_t* system, const uint count){
    printf("SOLUTION = ");
    for(uint i=0; i<count; i++){
        printf("%.6e ; ", solution[i]);
    }
    printf("\nSYSTEM VALUE = ");
    for(uint i=0; i<count; i++){
        printf("%.6e ; ", system[i](solution, count));
    }
    printf("\n");
}


void calc_shifted_value(const double shift, const double* delta,
                        const double* ini_vec, double* shifted_val, const uint count){
    memcpy(shifted_val, ini_vec, count*sizeof(double));
    cblas_daxpy(count, shift, delta, 1, shifted_val, 1);
}

double find_min_parameter(const double l_bnd, const double r_bnd, const double* delta,
                          const double* ini_vec, const function_t* system,
                          const uint count, const double param_precision){
    double cur_r_bnd = r_bnd;
    double cur_l_bnd = l_bnd;
    double mid_left_val[count];
    double mid_left_rhs[count];
    double mid_right_val[count];
    double mid_right_rhs[count];
    while(cur_r_bnd - cur_l_bnd>param_precision){
        double step = (cur_r_bnd - cur_l_bnd)/3;
        double mid_left = cur_l_bnd + step;
        double mid_right = cur_r_bnd - step;
        calc_shifted_value(mid_left, delta, ini_vec, mid_left_val, count);
        fill_func_val(system, mid_left_val, count, mid_left_rhs);
        double mid_left_defect = cblas_dnrm2(count, mid_left_rhs, 1);
        calc_shifted_value(mid_right, delta, ini_vec, mid_right_val, count);
        fill_func_val(system, mid_right_val, count, mid_right_rhs);
        double mid_right_defect = cblas_dnrm2(count, mid_right_rhs, 1);
        if(mid_left_defect<mid_right_defect){
            cur_r_bnd = mid_right;
        } else {
            cur_l_bnd = mid_left;
        }
    }
    if(cur_r_bnd==r_bnd) return r_bnd;
    return cur_r_bnd;
}

void newton_method(double* solution,const function_t* system, const uint count,
                    const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num){
    uint iter = 0;
    double rhs[count];
    fill_func_val(system, solution, count, rhs);
    double defect = cblas_dnrm2(count, rhs, 1);
    while(defect>precision){
        double der_matrix[count*count];
        fill_derivative_matrix(system, solution, count, der_delta, der_matrix);
        inverse_square_matrix(der_matrix, count);
        double delta_x[count];
        cblas_dgemv(CblasRowMajor, CblasNoTrans, count, count, -1.0, der_matrix,
                    count, rhs, 1, 0.0, delta_x, 1);
        double alpha = find_min_parameter(0.0, 1.0, delta_x, solution, system,
                                          count, 1e-2);
        cblas_daxpy(count, alpha, delta_x, 1, solution, 1);
        fill_func_val(system, solution, count, rhs);
        defect = cblas_dnrm2(count, rhs, 1);
        iter++;
        if(iter>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = iter;
}


void fill_norm_column(double* norm_column, const double* der_matrix, const uint count){
    for(uint i=0; i<count; i++){
        norm_column[i] = cblas_dnrm2(count, der_matrix+i*count, 1);
    }
}


void normalize_vector(double* vec, const double* norm_column, const uint count){
    for(uint i=0; i<count; i++){
        vec[i] /= norm_column[i];
    }
}

double calculate_normalized_defect(const double* rhs, const double* der_matrix, const uint count){
    double norm_column[count];
    fill_norm_column(norm_column, der_matrix, count);
    double tmp_rhs[count];
    memcpy(tmp_rhs, rhs, count*sizeof(double));
    normalize_vector(tmp_rhs, norm_column, count);
    double defect = cblas_dnrm2(count, tmp_rhs, 1);
    return defect;
}

void newton_method_normalized(double* solution,const function_t* system, const uint count,
                    const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num){
    uint iter = 0;
    double rhs[count];
    fill_func_val(system, solution, count, rhs);
    double general_defect = cblas_dnrm2(count, rhs, 1);
    double der_matrix[count*count];
    double delta_x[count];
    double tmp_solution[count];
    while(general_defect>precision){
        fill_derivative_matrix(system, solution, count, der_delta, der_matrix);
        double cur_norm_defect = calculate_normalized_defect(rhs, der_matrix, count);
        inverse_square_matrix(der_matrix, count);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, count, count, -1.0, der_matrix,
                    count, rhs, 1, 0.0, delta_x, 1);
        double alpha = 1.0;
        double next_norm_defect = 0;//updated later
        do{
            alpha -= 1e-2;
            memcpy(tmp_solution, solution, count*sizeof(double));
            cblas_daxpy(count, alpha, delta_x, 1, tmp_solution, 1);
            fill_func_val(system, tmp_solution, count, rhs);
            fill_derivative_matrix(system, tmp_solution, count, der_delta, der_matrix);
            next_norm_defect = calculate_normalized_defect(rhs, der_matrix, count);
        }while(next_norm_defect>cur_norm_defect && alpha>1e-2);
        cblas_dswap(count, solution, 1, tmp_solution, 1);
        general_defect = cblas_dnrm2(count, rhs, 1);
        iter++;
        if(iter>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = iter;
}



void newton_method_classic(double* solution,const function_t* system,
       const uint count, const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num){
    uint iter = 0;
    double rhs[count];
    fill_func_val(system, solution, count, rhs);
    double defect = cblas_dnrm2(count, rhs, 1);
    double der_matrix[count*count];
    while(defect>precision){
        fill_derivative_matrix(system, solution, count, der_delta, der_matrix);
        inverse_square_matrix(der_matrix, count);
        double delta_x[count];
        cblas_dgemv(CblasRowMajor, CblasNoTrans, count, count, -1.0, der_matrix,
                    count, rhs, 1, 0.0, delta_x, 1);
        cblas_daxpy(count, 1.0, delta_x, 1, solution, 1);
        fill_func_val(system, solution, count, rhs);
        defect = cblas_dnrm2(count, rhs, 1);
        iter++;
        if(iter>=max_iter) break;
    }
    if(iter_num!=NULL) *iter_num = iter;
}


void error_hendler(const int info){
    if(info>0){
        fprintf(stderr, "I have exactly ZERO in %i position\n", info);
        exit(EXIT_FAILURE);
    }
    else if(info<0){
        fprintf(stderr, "The %i-th argument has illegal value\n", -info);
        exit(EXIT_FAILURE);
    }
}

