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
void displayExecutionTime(clock_t, clock_t, TimeFormat_t);
void useGaussMethod(double **, double *, double *, unsigned int);
void displayMatrix(double **, double *, unsigned int);
void displaySolution(double *, unsigned int);

int main()
{
    unsigned int N = 3; // Number of equations

    // Memory allocation
    
    double **A = (double **)malloc(N * sizeof(double *));
    for (unsigned int i = 0; i < N; i++)
    {
        A[i] = (double *)malloc(N * sizeof(double));
    }
    double *Y = (double *)malloc(N * sizeof(double));
    double *X = (double *)malloc(N * sizeof(double));

    // Filling arrays
    
    A[0][0] = 2.351386;
    A[0][1] = 4.721297;
    A[0][2] = 1.983437;
    A[1][0] = 5.534478;
    A[1][1] = 2.375764;
    A[1][2] = 1.392935;
    A[2][0] = 2.894297;
    A[2][1] = 3.154358;
    A[2][2] = 4.653876;

    Y[0] = 36.313643;
    Y[1] = 47.645311;
    Y[2] = 37.287654;

    // Calculate and display

    printf("Source matrix:\n");
    displayMatrix(A, Y, N);

    clock_t execBegin = clock();
    useGaussMethod(A, Y, X, N);
    clock_t execEnd = clock();

    printf("Solution:\n");
    displaySolution(X, N);

    // Execution time

    TimeFormat_t format = MKS;
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

        for (unsigned int i = k + 1; i < N; i++)
        {
            if (abs(A[i][k]) > max)
            {
                max = abs(A[i][k]);
                index = i;
            }
        }

        // Rearranging strings

        if (max < eps)
        {
            printf("There are no non-zero diagonal elements.\n");
        }

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
