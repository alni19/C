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
void normal_sum0_init(double* array, int N);/* initilizes a normal array with sum 0*/

#endif
