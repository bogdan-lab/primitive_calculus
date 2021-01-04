#include "mse_utils.h"






int main(){
    double* test_arr;
    uint count;
    read_column_from_file("test.txt", &test_arr, &count, 10);

    for(uint i=0; i<count; i++){
        printf("%.3e ; ", test_arr[i]);
    }
    printf("\n");



    return 0;
}
