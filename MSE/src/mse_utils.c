#include "mse_utils.h"



void read_column_from_file(const char *filename, double** arr, uint* count,
                           const uint init_alloc){
    FILE* f_in = fopen(filename, "r");
    if(!f_in){
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    uint idx = 0;
    uint multiplier = 1;
    double* tmp_arr = calloc(init_alloc, sizeof(*tmp_arr));
    while(fscanf(f_in, "%lf", &tmp_arr[idx])!=EOF){
        idx++;
        if(idx==init_alloc*multiplier){
            multiplier *= 2;
            tmp_arr = realloc(tmp_arr, init_alloc*multiplier*sizeof(*tmp_arr));
        }
    }
    tmp_arr = realloc(tmp_arr, idx*sizeof(*tmp_arr));
    *arr = tmp_arr;
    *count = idx;
    fclose(f_in);
}


void calculate_combination_matrix(double* arr, const double* x, const uint comb_num,
                                  const uint count_x, const function_t comb_func){
    double comb_for_given_x[comb_num];
    for(uint i=0; i<count_x; i++){
        comb_func(comb_for_given_x, comb_num, x[i]);
        memcpy(arr+i*comb_num, comb_for_given_x, comb_num*sizeof(double));
    }
}

void polynomial_combination(double* res_arr, const uint comb_num, const double x){
    for(uint i=0; i<comb_num; i++){
        res_arr[i] = pow(x, i);
    }
}

void legandre_combination(double* res_arr, const uint comb_num, const double x){
    for(uint i=0; i<comb_num; i++){
        if(i==0){
            res_arr[i] = 1;
            continue;
        }
        if(i==1){
            res_arr[i] = x;
            continue;
        }
        res_arr[i] = 2.0*i/(i+1)*x*res_arr[i-1] - 1.0*i/(i+1)*res_arr[i-2];
    }
}


void regularize_lhs_matrix(double* res_arr, const uint comb_num, const uint count_x,
                           const double* lhs_A){
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, comb_num, comb_num, count_x,
                1.0, lhs_A, comb_num, lhs_A, comb_num, 0.0, res_arr, comb_num);
}

void regularize_rhs_vector(double* res_arr, const uint comb_num, const uint count_x,
                           const double* lhs_A, const double* rhs){
    cblas_dgemv(CblasRowMajor, CblasTrans, count_x, comb_num, 1.0, lhs_A, comb_num,
                rhs, 1, 0.0, res_arr, 1);
}

void check_for_diagonal_max(const double* arr, const uint count){
    printf("MAX ELEMENTS ARE IN : ");
    for(uint i=0; i<count; i++){
        uint max_c = 0;
        double max = fabs(arr[i*count + 0]);
        for(uint j=0; j<count; j++){
            double el =fabs(arr[i*count + j]);
            if(max<el){
                max_c = j;
                max = el;
            }
        }
        printf(" ( %u ; %u ) ", i, max_c);
    }
    printf("\n");
}


void print_matrix(const double* arr, const uint rc, const uint cc){
    uint total=cc*rc;
    for(uint i=0; i<total; i++){
        if(i%cc==0) printf("\n");
        printf("  %.4e  ", arr[i]);
    }
    printf("\n");
}

