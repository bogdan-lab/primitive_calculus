#include "mse_utils.h"






int main(){

    const uint comb_num = 3;

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

    double lhs_A[count_x*comb_num];
    calculate_combination_matrix(lhs_A, x, comb_num, count_x,
                                 //legandre_combination);
                                 polynomial_combination);

    double reg_lhs[comb_num*comb_num];
    regularize_lhs_matrix(reg_lhs, comb_num, count_x, lhs_A);
    print_matrix(reg_lhs, comb_num, comb_num);
    check_for_diagonal_max(reg_lhs, comb_num);
    double reg_rhs[comb_num];
    regularize_rhs_vector(reg_rhs, comb_num, count_x, lhs_A, y);


    return 0;
}
