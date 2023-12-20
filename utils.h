#ifndef _UTILS_H_
#define _UTILS_H_

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <string.h>

typedef enum 
{
    SEC, // Seconds
    MS,  // Milliseconds
    MKS, // Microseconds
} TimeFormat_t;

double drand(double, double);
void displayMatrix(double **, double *, unsigned int);
void displaySolution(double *, unsigned int);
void swap(double *, double *);

#endif
