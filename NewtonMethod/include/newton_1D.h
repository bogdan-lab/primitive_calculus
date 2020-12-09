#ifndef NEWTON_1D_H
#define NEWTON_1D_H

#include <stdlib.h>

//#define VERBOSE

typedef enum{
    LEFT_DERIVATIVE,
    MID_DERIVATIVE,
    RIGHT_DERIVATIVE
} DerivType;

double get_derivative(double (*func)(), const double x0, const double delta,
                      const DerivType der_type);

double pure_calc_newton_1D(double (*func)(double), const double x0, const double precision,
    const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_made);

double pure_anal_newton_1D(double (*func)(double), const double x0, const double precision,
                       double (*deriv)(double), const uint max_iter, uint* iter_made);

double modified_calc_newton_1D(double (*func)(double), const double x0, const double precision,
                       const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_made);

double modified_anal_newton_1D(double (*func)(double), const double x0, const double precision,
                       double (*deriv)(double), const uint max_iter, uint* iter_made);

double bisection(double (*func)(double), const double x_left, const double x_right,
                 const double precision, const uint max_iter, uint* iter_num);

#endif  //NEWTON_1D_H
