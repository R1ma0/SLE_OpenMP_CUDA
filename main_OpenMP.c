#include <omp.h>
#include "utils.h"

void useGaussMethod(double **, double *, double *, unsigned int);

int main()
{
    int maxThreads = omp_get_max_threads();
    unsigned int N = 1000; // Number of equations
                        
    printf("Max threads = %d\n", maxThreads);
    printf("Number of equations = %d\n", N);

    // Memory allocation
    
    unsigned int i;
    
    double **A = (double **)malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
    {
        A[i] = (double *)malloc(N * sizeof(double));
    }
    double *Y = (double *)malloc(N * sizeof(double));
    double *X = (double *)malloc(N * sizeof(double));

    // Filling arrays
    
    for (i = 0; i < N; i++)
    {
        for (unsigned int j = 0; j < N; j++)
        {
            A[i][j] = drand(0.0, 10.0);
        }

        Y[i] = drand(25.0, 56.0);
    }

    // Calculate and display

    //printf("# Source matrix:\n");
    //displayMatrix(A, Y, N);

    clock_t execBegin = clock();
    useGaussMethod(A, Y, X, N);
    clock_t execEnd = clock();

    //printf("# Solution:\n");
    //displaySolution(X, N);

    // Execution time
   
    TimeFormat_t format = MS;
    displayExecutionTime(execBegin, execEnd, format);

    // Freeing memory

    for (unsigned int i = 0; i < N; i++)
    {
        free(A[i]);
    }
    free(A);
    free(Y);
    free(X);

    return 0;
}

void useGaussMethod(double **A, double *Y, double *X, unsigned int N)
{
    const double eps = DBL_EPSILON;
    double max;
    unsigned int k = 0;
    unsigned int index;

    while (k < N)
    {
        // Search max A[i][k]

        max = abs(A[k][k]);
        index = k;
    
        unsigned int i;

#pragma omp parallel for private(i)
        for (i = k + 1; i < N; i++)
        {
            if (abs(A[i][k]) > max)
            {
                max = abs(A[i][k]);
                index = i;
            }
        }

        // Rearranging strings

        unsigned int j;

#pragma omp parallel for private(j)
        for (j = 0; j < N; j++)
        {
            double tmp = A[k][j];
            A[k][j] = A[index][j];
            A[index][j] = tmp;
        }

        double tmp = Y[k];
        Y[k] = Y[index];
        Y[index] = tmp;

        // Normalization

#pragma omp parallel for private(i, j, tmp)
        for (i = k; i < N; i++)
        {
            tmp = A[i][k];
            
            if (abs(tmp) < eps)
            {
                continue;
            }

            for (j = k; j < N; j++)
            {
                A[i][j] = A[i][j] / tmp;
            }
            
            Y[i] = Y[i] / tmp;

            if (i == k)
            {
                continue;
            }

            for (j = k; j < N; j++)
            {
                A[i][j] = A[i][j] - A[k][j];
            }

            Y[i] = Y[i] - Y[k];
        }

        k++;
    }

    // Reverse substitution

    int i,
        t;

    for (t = N - 1; t >= 0; t--)
    {
        X[t] = Y[t];

#pragma omp parallel for private(i)
        for (i = 0; i < t; i++)
        {
            Y[i] -= A[i][t] * X[t];
        }
    }
}
