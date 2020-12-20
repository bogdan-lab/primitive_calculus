#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "stdlib.h"
#include "lapacke.h"
#include "gsl/gsl_cblas.h"


typedef double(*function_t)(const double*, const uint);

void newton_method(double* solution, function_t* system, const uint count,
                    const double precision, const double der_delta,
                   const uint max_iter, uint* iter_num);


void inverse_square_matrix(double* matrix, const uint n);
void fill_func_val(function_t* system, const double* point, const uint count,
                   double* val_vec);
int check_solution(function_t* system, const double* solution, const uint count,
                   const double precision);
void fill_derivative_matrix(function_t* system, const double* point, const uint count,
                            const double der_delta, double* der_matrix);
double get_derivative_val(function_t func, const uint pos, const double* point,
                          const uint count, const double der_delta);

void error_hendler(const int info);
