#ifndef NEWTON_1D_H
#define NEWTON_1D_H


typedef enum{
    LEFT_DERIVATIVE,
    MID_DERIVATIVE,
    RIGHT_DERIVATIVE
} DerivType;

double solve_newton_1D(double (*func)(double), const double x0, const double precision,
    const DerivType der_type, const double der_delta, const uint max_iter, int* status);


#endif  //NEWTON_1D_H
