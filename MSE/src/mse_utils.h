#ifndef MSE_UTILS_H
#define MSE_UTILS_H



#include<stdlib.h>
#include<stdio.h>


void read_column_from_file(const char* filename, double** arr, uint* count,
                           const uint init_alloc);

#endif //MSE_UTILS_H
