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

double get_second_derivative(double (*func)(double), const double x,
                             const double der_delta);

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

double secant_method(double (*func)(double), const double x_0, const double x_1,
                 const double precision, const uint max_iter, uint* iter_num);

double chord_secant_method(double (*func)(double), const double x_0, const double x_1,
                       const double precision, const uint max_iter, uint* iter_num);

double modified_secant_method(double (*func)(double), const double x_0, const double x_1,
                 const double precision, const uint max_iter, uint* iter_num);

void swap(double* lhs, double* rhs);

double third_order_newton(double (*func)(double), const double x0, const double precision,
  const DerivType der_type, const double der_delta, const uint max_iter, uint* iter_num);

double third_order_newton_anal(double (*func)(double), const double x0,
  const double precision, double (*deriv)(double), double (*deriv_sec)(double),
       const uint max_iter, uint* iter_num);

double third_order_newton_modified(double (*func)(double), const double x0,
        const double precision, const DerivType der_type, const double der_delta,
           const uint max_iter, uint* iter_num);

double third_order_newton_modified_anal(double (*func)(double), const double x0,
     const double precision, double (*deriv)(double), double (*deriv_sec)(double),
           const uint max_iter, uint* iter_num);

#endif  //NEWTON_1D_H
