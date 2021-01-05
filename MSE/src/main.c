#include "mse_utils.h"


#include<lapacke.h>

void run_mse(const double* x, const uint count_x, const double* y,
             const uint comb_num, const function_t comb_func, const char* fout_name){
    double lhs_A[count_x*comb_num];
    calculate_combination_matrix(lhs_A, x, comb_num, count_x,
                                 comb_func);

    double reg_lhs[comb_num*comb_num];
    regularize_lhs_matrix(reg_lhs, comb_num, count_x, lhs_A);
    //print_matrix(reg_lhs, comb_num, comb_num);
    check_for_diagonal_max(reg_lhs, comb_num);
    double reg_rhs[comb_num];
    regularize_rhs_vector(reg_rhs, comb_num, count_x, lhs_A, y);

    double weights[comb_num];
    solve_equation(weights, reg_lhs, reg_rhs, comb_num);
    save_weight_to_file(fout_name, weights, comb_num);
}

int main(){

    const uint comb_num = 30;

    double* x;
    uint count_x;
    double* y;
    uint count_y;

    read_column_from_file("x.txt", &x, &count_x, 10);
    read_column_from_file("y.txt", &y, &count_y, 10);
    cblas_dnrm2(count_x, x, 1);
    if(count_x!=count_y){
        fprintf(stderr, "Incorrect input, count_x (%u) != count_y (%u)",
                count_x, count_y);
        exit(EXIT_FAILURE);
    }

    run_mse(x, count_x, y, comb_num, polynomial_combination, "weight_pol.txt");
    run_mse(x, count_x, y, comb_num, legandre_combination, "weight_leg.txt");

    return 0;
}
