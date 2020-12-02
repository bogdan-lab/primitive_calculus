#include <stdlib.h>
#include <stdio.h>
#include "newton_1D.h"

double solve_newton_1D(double (*func)(double), const double x0, const double precision,
                       const DerivType der_type, const double der_delta, const uint max_iter, int* status){
    printf("**** solve_newton_1D ****\n");
    double x = x0;
    double defect = func(x);
    double xmin = x;
    double defec_min = defect;
    uint idx = 0;
    while(abs(defect)>precision){
        double derivative = 0;
        switch (der_type) {
        case LEFT_DERIVATIVE:
            derivative = (defect - func(x-der_delta))/der_delta;
            break;
        case RIGHT_DERIVATIVE:
            derivative = (func(x+der_delta)-defect)/der_delta;
            break;
        case MID_DERIVATIVE:;
            derivative = (func(x+der_delta) - func(x-der_delta))/der_delta;
            break;
        default:
            fprintf(stderr, "Incorrect derivative type flag\n");
            exit(EXIT_FAILURE);
        }
        x -= defect/derivative;
        idx++;
        printf("iteration[%u]\tdelta = %.6e\n", idx, defect);
        xmin = defec_min<defect ? x : xmin;
        if(idx>=max_iter){
            printf("Limit of the iteration number was reached!\n");
            *status = -1;
            return xmin;
        }
    }
    *status = 1;
    return x;
}

