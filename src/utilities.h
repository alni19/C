#ifndef UTILITIES

#define UTILITIES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NOMEMORY 100

char does_it_exist(char *filename);
void gotosleep(int numseconds);
void erasefile(char *filename);
double drawnormal(void);
void split(double* price,char* str);
void normal_sum0_init(double* arr, int N);/* initilizes a normal array with sum 0*/
void calculate_v(const double** matrix, double* v,int numassets, int numdays);
void perturb(const double** matrix, double** target, int numassets, int numdays,const double* v);

#endif
