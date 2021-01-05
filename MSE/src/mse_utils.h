#ifndef MSE_UTILS_H
#define MSE_UTILS_H



#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<gsl/gsl_cblas.h>
#include<lapacke.h>

typedef void(*function_t)(double* res_arr, const uint comb_count, const double x);

void read_column_from_file(const char* filename, double** arr, uint* count,
                           const uint init_alloc);
void calculate_combination_matrix(double* arr, const double* x, const uint comb_num,
                                  const uint count_x, const function_t comb_func);

void polynomial_combination(double* res_arr, const uint comb_num, const double x);

void legandre_combination(double* res_arr, const uint comb_num, const double x);
void regularize_lhs_matrix(double* res_arr, const uint comb_num, const uint count_x,
                           const double* lhs_A);
void regularize_rhs_vector(double* res_arr, const uint comb_num, const uint count_x,
                           const double* lhs_A, const double* rhs);
void print_matrix(const double* arr, const uint rc, const uint cc);
void check_for_diagonal_max(const double* arr, const uint count);
void solve_equation(double* x, const double* A, const double* b, const uint count);
void error_hendler(const int info);
void inverse_square_matrix(double* matrix, const uint n);
void save_weight_to_file(const char* filename, const double* weights, const uint count);
#endif //MSE_UTILS_H
