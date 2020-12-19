#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "lapacke.h"


void error_hendler(const int info);
void inverse_square_matrix(double* matrix, const uint n);
void matrix_vector_mult(double** solution, const double* inv_m, const double* rhs,
                   const uint n);
void randomly_fill_array(double** arr, const uint size);
void check_solution(const double* solution, const double* lhs, const double* rhs,
                    const uint eq_num, const double defect);

void my_inverse_2x2(double** result, const double* input);

int main(){
    srand((uint)time(NULL));
    const double delta = 1e-10;
    const uint eq_num = 10;
    const uint matrix_size = eq_num*eq_num;
    double* rhs;
    double* lhs;
    randomly_fill_array(&rhs, eq_num);
    randomly_fill_array(&lhs, matrix_size);
    double* lhs_inv = calloc(matrix_size, sizeof(*lhs_inv));
    memcpy(lhs_inv, lhs, matrix_size*sizeof(double));
    inverse_square_matrix(lhs_inv, eq_num);
    double* solution;
    matrix_vector_mult(&solution, lhs_inv, rhs, eq_num);
    check_solution(solution, lhs, rhs, eq_num, delta);
    return 0;
}


void error_hendler(const int info){
    if(info!=0){
        fprintf(stderr, "Something went wrong!");
        exit(EXIT_FAILURE);
    }
}


void inverse_square_matrix(double* matrix, const uint n){
    int ipiv[n];
    int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, matrix, n, ipiv);
    error_hendler(info);
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, matrix, n, ipiv);
    error_hendler(info);
}


void matrix_vector_mult(double** solution, const double* matrix, const double* vec,
                   const uint n){
    double* res = calloc(n, sizeof(*res));
    for(uint i=0; i<n; i++){
        res[i] = 0;
        for(uint j=0; j<n; j++){
            res[i] += matrix[n*i+j]*vec[j];
        }
    }
    *solution = res;
}


void randomly_fill_array(double** arr, const uint size){
    double* tmp_arr = calloc(size, sizeof(*tmp_arr));
    for(uint i=0; i<size; i++){
        tmp_arr[i] = 1.0*rand()/RAND_MAX;
    }
    *arr = tmp_arr;
}

void check_solution(const double* solution, const double* lhs, const double* rhs,
                    const uint eq_num, const double delta){
    double* res;
    matrix_vector_mult(&res, lhs, solution, eq_num);
    for(uint i=0; i<eq_num; i++){
        double defect = fabs(rhs[i] - res[i]);
        if(defect>delta){
            fprintf(stderr, "Row %i has incorrect result : defect = %.6e\n", i, defect);
            exit(EXIT_FAILURE);
        }
    }
    printf("Check was succesfull!\n");
}

void my_inverse_2x2(double** result, const double* input){
    double det = input[0]*input[3]-input[1]*input[2];
    double* tmp = calloc(4, sizeof(*tmp));
    tmp[0] = input[3]/det;
    tmp[1] = -input[1]/det;
    tmp[2] = -input[2]/det;
    tmp[3] = input[0]/det;
    *result = tmp;
}
