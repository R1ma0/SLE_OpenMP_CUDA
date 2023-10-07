#ifndef _UTILS_H_
#define _UTILS_H_

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

double calcCPUTimeUsage(clock_t, clock_t, double);
double drand(double, double);
void displayExecutionTime(clock_t, clock_t, TimeFormat_t);
void displayMatrix(double **, double *, unsigned int);
void displaySolution(double *, unsigned int);

#endif
