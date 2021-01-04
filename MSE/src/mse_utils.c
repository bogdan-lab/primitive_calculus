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
