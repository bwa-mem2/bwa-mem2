#include "sorted_doubles_rmi.h"
#include <math.h>
#include <cmath>
#include <fstream>
//#include <filesystem>
#include <iostream>
#include <stdio.h>
#include <sys/types.h>
#include<sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[])
{
    char *filename = argv[1];
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        printf("%s file not found\n", filename);
        exit(0);
    }
    printf("query file good\n");
    int64_t n = 0;
    infile.read((char*)&n, sizeof(int64_t));
    printf("n = %ld\n", n);
    double *tmp_data = (double *)malloc(n * sizeof(double));
    if (tmp_data == NULL)
    {
        printf("tmp_data is NULL\n");
        exit(0);
    }
    infile.read((char*)tmp_data, n * sizeof(double));
    printf("read query file\n");

    int64_t i;
    sorted_doubles_rmi::load("rmi_data");
    printf("loaded rmi data\n");
    sorted_doubles_rmi::save(argv[2]);
    double total_log_err = 0;
    double total_log_gap = 0;
    int64_t max_err = 0;
    int64_t max_gap = 0;
    int64_t err_hist[20];
    memset(err_hist, 0, 20*8);
    for(i = 0; i < n; i++)
    {
        size_t err = 0;
        size_t pos = sorted_doubles_rmi::lookup(tmp_data[i], &err);
        //printf("%d] key = %E, %lu, %lu, %lu\n", i, tmp_data[i], pos, err, labs(pos-i));
        size_t gap = labs(pos-i);
        if(gap > err)
            printf("%d] key = %E, %lu, %lu, %lu\n", i, tmp_data[i], pos, err, labs(pos-i));

        int64_t log2_err = log2(1.0 * err);
        err_hist[log2_err]++;
        if(max_err < err) max_err = err;
        if(max_gap < gap) max_gap = gap;
        total_log_err += log2(1.0 + err);
        total_log_gap += log2(1.0 + gap);
    }
    printf("avg_log2_err = %lf, avg_log2_gap = %lf\n", total_log_err / n, total_log_gap / n);
    printf("max_err = %ld, max_gap = %ld\n", max_err, max_gap);
    for(i = 0; i < 20; i++)
    {
        printf("%ld] log2_err freq = %ld\n", i, err_hist[i]);
    }
    sorted_doubles_rmi::cleanup();
}
