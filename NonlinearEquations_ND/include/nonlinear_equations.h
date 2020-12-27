#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "stdlib.h"
#include "lapacke.h"
#include "gsl/gsl_cblas.h"


typedef double(*function_t)(const double*, const uint);

void newton_method(double* solution, const function_t* system, const uint count,
                    const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num);

void newton_method_classic(double* solution,const function_t* system,
       const uint count, const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num);

void newton_method_normalized(double* solution,const function_t* system, const uint count,
                    const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num);

void inverse_square_matrix(double* matrix, const uint n);
void fill_func_val(const function_t* system, const double* point, const uint count,
                   double* val_vec);
int check_solution(const function_t* system, const double* solution, const uint count,
                   const double precision);
void fill_derivative_matrix(const function_t* system, const double* point, const uint count,
                            const double der_delta, double* der_matrix);
double get_derivative_val(const function_t func, const uint pos, const double* point,
                          const uint count, const double der_delta);
void print_result(const double* solution, const function_t* system, const uint count);
void error_hendler(const int info);

void normalize_vector(double* vec, const double* norm_column, const uint count);
void fill_norm_column(double* norm_column, const double* der_matrix, const uint count);

double calculate_normalized_defect(const double* rhs, const double* der_matrix, const uint count);

double find_min_parameter(const double l_bnd, const double r_bnd, const double* delta,
                          const double* ini_vec, const function_t* system,
                          const uint count, const double param_precision);
void calc_shifted_value(const double shift, const double* delta,
                        const double* ini_vec, double* shifted_val, const uint count);
