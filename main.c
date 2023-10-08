#include "utils.h"

void useGaussMethod(double **, double *, double *, unsigned int);

int main(int argc, char *argv[])
{
    unsigned int N; // Number of equations

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

    unsigned int j;

#pragma omp parallel for schedule(dynamic, 32) private(j)
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i][j] = j;
        }

        Y[i] = i;
    }

    // Calculate and display

    double execBegin = omp_get_wtime();
    useGaussMethod(A, Y, X, N);
    double execEnd = omp_get_wtime();

    // Execution time

    printf("Execution time: %f sec\n", execEnd - execBegin);

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

        for (unsigned int i = k + 1; i < N; i++)
        {
            if (abs(A[i][k]) > max)
            {
                max = abs(A[i][k]);
                index = i;
            }
        }

        // Rearranging strings

        for (unsigned int j = 0; j < N; j++)
        {
            double tmp = A[k][j];
            A[k][j] = A[index][j];
            A[index][j] = tmp;
        }

        double tmp = Y[k];
        Y[k] = Y[index];
        Y[index] = tmp;

        // Normalization

        for (unsigned int i = k; i < N; i++)
        {
            double tmp = A[i][k];
            
            if (abs(tmp) < eps)
            {
                continue;
            }

            for (unsigned int j = k; j < N; j++)
            {
                A[i][j] = A[i][j] / tmp;
            }
            
            Y[i] = Y[i] / tmp;

            if (i == k)
            {
                continue;
            }

            for (unsigned int j = k; j < N; j++)
            {
                A[i][j] = A[i][j] - A[k][j];
            }

            Y[i] = Y[i] - Y[k];
        }

        k++;
    }

    // Reverse substitution

    for (int t = N - 1; t >= 0; t--)
    {
        X[t] = Y[t];
        
        for (int i = 0; i < t; i++)
        {
            Y[i] -= A[i][t] * X[t];
        }
    }
}
