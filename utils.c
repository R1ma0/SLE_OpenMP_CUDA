#include "utils.h"

double drand(double min, double max)
{
    return ((double)rand() * (max - min)) / (double)RAND_MAX + min;
}

void displayMatrix(double **A, double *Y, unsigned int N)
{
    for (unsigned int i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            printf("%f*x", A[i][j]);

            if (j < N - 1)
            {
                printf(" + ");
            }
        }

        printf(" = %f\n", Y[i]);
    }
}

void displaySolution(double *X, unsigned int N)
{
    for (unsigned int i = 0; i < N; i++)
    {
        printf("X%d = %f\n", i + 1, X[i]);
    }
}

double calcCPUTimeUsage(clock_t begin, clock_t end, double divider)
{
    return ((double)(end - begin)) / (CLOCKS_PER_SEC / divider);
}

void displayExecutionTime(clock_t begin, clock_t end, TimeFormat_t format)
{
    double divider;
    char title[4];

    switch (format)
    {
        case SEC:
            divider = 1.0;
            strcpy(title, "sec");
            break;
        case MS:
            divider = 1000.0;
            strcpy(title, "ms");
            break;
        case MKS:
            divider = 1000000.0;
            strcpy(title, "mks");
            break;
        default:
            divider = 1.0;
            strcpy(title, "sec");
    }

    double time = calcCPUTimeUsage(begin, end, divider);
    printf("Execution time: %f %s\n", time, title);
}
