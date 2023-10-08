#include <omp.h>
#include "utils.h"

void useGaussMethod(double **, double *, double *, unsigned int);

int main(int argc, char *argv[])
{
    unsigned int N; // Number of equations
     
    printf("=== Parallel algorithm ===\n");

    if (argc <= 1)
    {
        printf("Enter size of the matrix.\n");
        return 0;
    }
    else
    {
        N = atoi(argv[1]); 
        printf("Number of equations = %d\n", N);
    }
                        
    int maxThreads = omp_get_max_threads();
    printf("Max threads = %d\n", maxThreads);

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

    // Calculate 

    clock_t execBegin = clock();
    useGaussMethod(A, Y, X, N);
    clock_t execEnd = clock();

    // Execution time
   
    TimeFormat_t format = SEC;
    displayExecutionTime(execBegin, execEnd, format);

    // Checking the correctness of the solution
    
    unsigned int M = 3;
    unsigned int k;
    
    double *correctSolution = (double *)malloc(M * sizeof(double));
    correctSolution[0] = 7.0;
    correctSolution[1] = 5.0;
    correctSolution[2] = 2.0;

    double **a = (double **)malloc(M * sizeof(double *));
    for (k = 0; k < M; k++)
    {
        a[k] = (double *)malloc(M * sizeof(double));
    }
    double *y = (double *)malloc(M * sizeof(double));
    double *x = (double *)malloc(M * sizeof(double));

    a[0][0] = 2.0;
    a[0][1] = 4.0;
    a[0][2] = 1.0;
    a[1][0] = 5.0;
    a[1][1] = 2.0;
    a[1][2] = 1.0;
    a[2][0] = 2.0;
    a[2][1] = 3.0;
    a[2][2] = 4.0;

    y[0] = 36.0;
    y[1] = 47.0;
    y[2] = 37.0;

    printf("\n=== Checking the correctness of the algorithm  ===\n");

    printf("Test matrix:\n");
    displayMatrix(a, y, M);

    useGaussMethod(a, y, x, M);

    printf("Resultion solution:\n");
    displaySolution(x, M);

    printf("Correct solution\n");
    displaySolution(correctSolution, M);

    // Freeing memory

    for (i = 0; i < N; i++)
    {
        free(A[i]);
    }
    for (k = 0; k < M; k++)
    {
        free(a[k]);
    }
    free(A);
    free(a);
    free(Y);
    free(y);
    free(X);
    free(x);
    free(correctSolution);

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

#pragma omp parallel for schedule(static) private(i)
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
        double tmp;

#pragma omp parallel for schedule(static) private(j, tmp)
        for (j = 0; j < N; j++)
        {
            tmp = A[k][j];
            A[k][j] = A[index][j];
            A[index][j] = tmp;
        }

        tmp = Y[k];
        Y[k] = Y[index];
        Y[index] = tmp;

        // Normalization

#pragma omp parallel for schedule(static) private(i, j)
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

#pragma omp parallel for schedule(static) private(t)
    for (t = N - 1; t >= 0; t--)
    {
        X[t] = Y[t];

        for (i = 0; i < t; i++)
        {
            Y[i] -= A[i][t] * X[t];
        }
    }
}
