#include "utils.h"

double drand(double min, double max)
{
    srand(15);
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
