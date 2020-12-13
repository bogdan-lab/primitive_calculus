#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "nonlinear_equation_1D.h"

const double precision = 1e-9;
const uint max_iter = 100;
const double der_delta = 1e-6;
static uint f_call_count = 0;


typedef struct Result{
    double all_time;
    double avg_func_calls;
    uint min_func_calls;
    uint max_func_calls;
    double avg_iter_num;
    uint min_iter_num;
    uint max_iter_num;
    uint successful_tests;
    uint nan_failed_tests;
    uint iter_limited_tests;
}Result;

void print_result(const char* method_name, const Result* res);

void StartTimer(time_t* start);
double GetDuration(time_t* start, time_t* end);

double get_avg(const uint* arr, const uint count);
double get_min(const uint* arr, const uint count);
double get_max(const uint* arr, const uint count);

double func_1(const double x);
double func_1_deriv(const double x);
double func_1_deriv2(const double x);

double func_2(const double x);
double func_2_deriv(const double x);
double func_2_deriv2(const double x);


double func_3(const double x);
double func_3_deriv(const double x);
double func_3_deriv2(const double x);

double func_4(const double);

Result benchmark_anal_newton(const uint test_num, const double l_guess_bnd,
      const double r_guess_bnd, double(*func)(double), double(*deriv)(double));

Result benchmark_calc_newton(const uint test_num, const double l_guess_bnd,
      const double r_guess_bnd, double(*func)(double));

Result benchmark_anal_newton_3ord(const uint test_num, const double l_guess_bnd,
      const double r_guess_bnd, double(*func)(double), double(*deriv)(double),
                                  double(*deriv2)(double));

Result benchmark_modif_secant(const uint test_num, const double l_guess_bnd,
      const double r_guess_bnd, double(*func)(double));



void test_func_1();
void test_func_2();
void test_func_3();
void test_func_4();


int main(){
    //test_func_1();
    //--> third order works badsince it is very sensitive for function (f^2, f'^{-3}) - make it unstable
    // newton and secant are close in iteration number, but secant requires less function calls
    //test_func_2();
    // --> third order is goood here - the best iteration score. Secant - best function calls!
    //test_func_3();
    //second order derivative is inf in solution, but other multipliers compensate it!
    //third order method is not that bad here, but classical newton is better.
    //secant is still close to newton and has less function calls, but slightly more iterations
    test_func_4();
    //As expected in case when derivative is crazy secant method is much better than newton
    //about 2 times less iterations and several times less function calls
    return 0;
}


void StartTimer(time_t* start){
    f_call_count = 0;
    *start = clock();
}


double GetDuration(time_t* start, time_t* end){
    *end = clock();
    double T = (double)(*end - *start)/(double)(CLOCKS_PER_SEC);
    return T;
}

double func_1(const double x){
    f_call_count++;
    return x+2-exp(x);
    //x1 = -1.841405660; x2 = 1.146193221
}

double func_1_deriv(const double x){
    f_call_count++;
    return 1-exp(x);
}

double func_1_deriv2(const double x){
    f_call_count++;
    return -exp(x);
}


double func_2(const double x){
    f_call_count++;
    return x*x;
}

double func_2_deriv(const double x){
    f_call_count++;
    return 2*x;
}

double func_2_deriv2(const double x){
    f_call_count++;
    isfinite(x);
    return 2;
}

double func_3(const double x){
    f_call_count++;
    return x + x*cbrt(x);
}

double func_3_deriv(const double x){
    f_call_count++;
    return 1 + 4.0/3.0*cbrt(x);
}

double func_3_deriv2(const double x){
    f_call_count++;
    return 4.0/9.0*1/cbrt(x*x);
}

double func_4(const double x){
    f_call_count++;
    if(x==0.0){return 0.0;}
    return x + x*x*sin(2.0/x);
}

Result benchmark_anal_newton(const uint test_num, const double l_guess_bnd,
   const double r_guess_bnd, double (*func)(double), double(*deriv)(double)){
    uint* iter_stat = calloc(test_num, sizeof(*iter_stat));
    uint* calls_stat = calloc(test_num, sizeof(*calls_stat));
    double time = 0;
    uint nan_failed_tests = 0;
    uint not_finished_tests = 0;
    time_t start, end;
    uint idx = 0;
    StartTimer(&start);
    for(uint i=0; i<test_num; i++){
        double x = l_guess_bnd + (r_guess_bnd - l_guess_bnd)*1.0*rand()/RAND_MAX;
        f_call_count = 0;
        double res = modified_anal_newton_1D(func, x, precision, deriv, max_iter,
                                         &iter_stat[idx]);
        calls_stat[idx] = f_call_count;
        if(!isfinite(res) || !isfinite(func(res))){
            nan_failed_tests++;
            continue;
        }
        if(iter_stat[idx]==max_iter) {
            not_finished_tests++;
            continue;
        }
        idx++;
    }
    time = GetDuration(&start, &end);
    iter_stat = realloc(iter_stat, idx*sizeof(*iter_stat));
    calls_stat = realloc(calls_stat, idx*sizeof(*calls_stat));
    Result stat;
    stat.all_time = time;
    stat.iter_limited_tests = not_finished_tests;
    stat.nan_failed_tests = nan_failed_tests;
    stat.successful_tests = test_num - nan_failed_tests - not_finished_tests;
    stat.avg_func_calls = get_avg(calls_stat, idx);
    stat.min_func_calls = get_min(calls_stat, idx);
    stat.max_func_calls = get_max(calls_stat, idx);
    stat.avg_iter_num = get_avg(iter_stat, idx);
    stat.min_iter_num = get_min(iter_stat, idx);
    stat.max_iter_num = get_max(iter_stat, idx);
    free(iter_stat);
    free(calls_stat);
    return stat;
}

Result benchmark_calc_newton(const uint test_num, const double l_guess_bnd,
   const double r_guess_bnd, double (*func)(double)){
    uint* iter_stat = calloc(test_num, sizeof(*iter_stat));
    uint* calls_stat = calloc(test_num, sizeof(*calls_stat));
    double time = 0;
    uint nan_failed_tests = 0;
    uint not_finished_tests = 0;
    time_t start, end;
    uint idx = 0;
    StartTimer(&start);
    for(uint i=0; i<test_num; i++){
        double x = l_guess_bnd + (r_guess_bnd - l_guess_bnd)*1.0*rand()/RAND_MAX;
        f_call_count = 0;
        double res = modified_calc_newton_1D(func, x, precision, MID_DERIVATIVE,
                   der_delta, max_iter, &iter_stat[idx]);
        calls_stat[idx] = f_call_count;
        if(!isfinite(res) || !isfinite(func(res))){
            nan_failed_tests++;
            continue;
        }
        if(iter_stat[idx]==max_iter) {
            not_finished_tests++;
            continue;
        }
        idx++;
    }
    time = GetDuration(&start, &end);
    iter_stat = realloc(iter_stat, idx*sizeof(*iter_stat));
    calls_stat = realloc(calls_stat, idx*sizeof(*calls_stat));
    Result stat;
    stat.all_time = time;
    stat.iter_limited_tests = not_finished_tests;
    stat.nan_failed_tests = nan_failed_tests;
    stat.successful_tests = test_num - nan_failed_tests - not_finished_tests;
    stat.avg_func_calls = get_avg(calls_stat, idx);
    stat.min_func_calls = get_min(calls_stat, idx);
    stat.max_func_calls = get_max(calls_stat, idx);
    stat.avg_iter_num = get_avg(iter_stat, idx);
    stat.min_iter_num = get_min(iter_stat, idx);
    stat.max_iter_num = get_max(iter_stat, idx);
    free(iter_stat);
    free(calls_stat);
    return stat;
}


Result benchmark_anal_newton_3ord(const uint test_num, const double l_guess_bnd,
   const double r_guess_bnd, double (*func)(double), double(*deriv)(double),
                                  double (*deriv2)(double)){
    uint* iter_stat = calloc(test_num, sizeof(*iter_stat));
    uint* calls_stat = calloc(test_num, sizeof(*calls_stat));
    double time = 0;
    uint nan_failed_tests = 0;
    uint not_finished_tests = 0;
    time_t start, end;
    uint idx = 0;
    StartTimer(&start);
    for(uint i=0; i<test_num; i++){
        double x = l_guess_bnd + (r_guess_bnd - l_guess_bnd)*1.0*rand()/RAND_MAX;
        f_call_count = 0;
        double res = third_order_newton_anal(func, x, precision, deriv,
                       deriv2, max_iter, &iter_stat[idx]);
        calls_stat[idx] = f_call_count;
        if(!isfinite(res) || !isfinite(func(res))){
            nan_failed_tests++;
            continue;
        }
        if(iter_stat[idx]==max_iter) {
            not_finished_tests++;
            continue;
        }
        idx++;
    }
    time = GetDuration(&start, &end);
    iter_stat = realloc(iter_stat, idx*sizeof(*iter_stat));
    calls_stat = realloc(calls_stat, idx*sizeof(*calls_stat));
    Result stat;
    stat.all_time = time;
    stat.iter_limited_tests = not_finished_tests;
    stat.nan_failed_tests = nan_failed_tests;
    stat.successful_tests = test_num - nan_failed_tests - not_finished_tests;
    stat.avg_func_calls = get_avg(calls_stat, idx);
    stat.min_func_calls = get_min(calls_stat, idx);
    stat.max_func_calls = get_max(calls_stat, idx);
    stat.avg_iter_num = get_avg(iter_stat, idx);
    stat.min_iter_num = get_min(iter_stat, idx);
    stat.max_iter_num = get_max(iter_stat, idx);
    free(iter_stat);
    free(calls_stat);
    return stat;
}



Result benchmark_modif_secant(const uint test_num, const double l_guess_bnd,
   const double r_guess_bnd, double (*func)(double)){
    uint* iter_stat = calloc(test_num, sizeof(*iter_stat));
    uint* calls_stat = calloc(test_num, sizeof(*calls_stat));
    double time = 0;
    uint nan_failed_tests = 0;
    uint not_finished_tests = 0;
    time_t start, end;
    uint idx = 0;
    StartTimer(&start);
    for(uint i=0; i<test_num; i++){
        double x0 = l_guess_bnd + (r_guess_bnd - l_guess_bnd)*1.0*rand()/RAND_MAX;
        double x1 = l_guess_bnd + (r_guess_bnd - l_guess_bnd)*1.0*rand()/RAND_MAX;
        f_call_count = 0;
        double res = modified_secant_method(func, x0, x1, precision, max_iter,
                                         &iter_stat[idx]);
        calls_stat[idx] = f_call_count;
        if(!isfinite(res) || !isfinite(func(res))){
            nan_failed_tests++;
            continue;
        }
        if(iter_stat[idx]==max_iter) {
            not_finished_tests++;
            continue;
        }
        idx++;
    }
    time = GetDuration(&start, &end);
    iter_stat = realloc(iter_stat, idx*sizeof(*iter_stat));
    calls_stat = realloc(calls_stat, idx*sizeof(*calls_stat));
    Result stat;
    stat.all_time = time;
    stat.iter_limited_tests = not_finished_tests;
    stat.nan_failed_tests = nan_failed_tests;
    stat.successful_tests = test_num - nan_failed_tests - not_finished_tests;
    stat.avg_func_calls = get_avg(calls_stat, idx);
    stat.min_func_calls = get_min(calls_stat, idx);
    stat.max_func_calls = get_max(calls_stat, idx);
    stat.avg_iter_num = get_avg(iter_stat, idx);
    stat.min_iter_num = get_min(iter_stat, idx);
    stat.max_iter_num = get_max(iter_stat, idx);
    free(iter_stat);
    free(calls_stat);
    return stat;
}



double get_avg(const uint* arr, const uint count){
    uint sum = 0;
    for(uint i=0; i<count; i++){
        sum += arr[i];
    }
    return 1.0*sum/count;
}

double get_min(const uint* arr, const uint count){
    if(count==0){
        fprintf(stderr, "Searching for min in empty array!\n");
        exit(EXIT_FAILURE);
    }
    uint min_val = arr[0];
    for(uint i=0; i<count; i++){
        min_val = fmin(min_val, arr[i]);
    }
    return min_val;
}

double get_max(const uint* arr, const uint count){
    if(count==0){
        fprintf(stderr, "Searching for max in empty array!\n");
        exit(EXIT_FAILURE);
    }
    uint max_val = arr[0];
    for(uint i=0; i<count; i++){
        max_val = fmax(max_val, arr[i]);
    }
    return max_val;
}

void print_result(const char* method_name, const Result* res){
    printf("**** %s ****\n", method_name);
    printf("FULL TIME [s] = %.6e\n", res->all_time);
    printf("TEST STAT: \n\tSucces = %u;\n\tNot finished = %u\n\tNaN failed = %u\n",
           res->successful_tests, res->iter_limited_tests, res->nan_failed_tests);
    printf("FUNC CALLS: \n\tmin = %u\n\tavg = %.2lf\n\tmax = %u\n",
           res->min_func_calls, res->avg_func_calls, res->max_func_calls);
    printf("ITERATIONS: \n\tmin = %u\n\tavg = %.2lf\n\tmax = %u\n",
           res->min_iter_num, res->avg_iter_num, res->max_iter_num);
    printf("==================================================================\n");
}


void test_func_1(){
    const uint tst_num = 1000000;
    const double xl = -10.0;
    const double xr = 10.0;
    Result res1 = benchmark_anal_newton(tst_num, xl, xr, func_1, func_1_deriv);
    print_result("Analitycal modified newton SECOND order", &res1);
    Result res2 = benchmark_anal_newton_3ord(tst_num, xl, xr, func_1,
                                             func_1_deriv, func_1_deriv2);
    print_result("Analitycal modified newton THIRD order", &res2);
    Result res3 = benchmark_modif_secant(tst_num, xl, xr, func_1);
    print_result("Modified Secant method", &res3);
}

void test_func_2(){
    const uint tst_num = 1000000;
    const double xl = -10.0;
    const double xr = 10.0;
    Result res1 = benchmark_anal_newton(tst_num, xl, xr, func_2, func_2_deriv);
    print_result("Analitycal modified newton SECOND order", &res1);
    Result res2 = benchmark_anal_newton_3ord(tst_num, xl, xr, func_2,
                                             func_2_deriv, func_2_deriv2);
    print_result("Analitycal modified newton THIRD order", &res2);
    Result res3 = benchmark_modif_secant(tst_num, xl, xr, func_2);
    print_result("Modified Secant method", &res3);
}

void test_func_3(){
    const uint tst_num = 1000000;
    const double xl = -10.0;
    const double xr = 10.0;
    Result res1 = benchmark_anal_newton(tst_num, xl, xr, func_3, func_3_deriv);
    print_result("Analitycal modified newton SECOND order", &res1);
    Result res2 = benchmark_anal_newton_3ord(tst_num, xl, xr, func_3,
                                             func_3_deriv, func_3_deriv2);
    print_result("Analitycal modified newton THIRD order", &res2);
    Result res3 = benchmark_modif_secant(tst_num, xl, xr, func_3);
    print_result("Modified Secant method", &res3);
}

void test_func_4(){
    const uint tst_num = 1000000;
    const double xl = -10.0;
    const double xr = 10.0;
    Result res1 = benchmark_calc_newton(tst_num, xl, xr, func_4);
    print_result("Calculus modified newton SECOND order", &res1);
    Result res3 = benchmark_modif_secant(tst_num, xl, xr, func_4);
    print_result("Modified Secant method", &res3);
}
