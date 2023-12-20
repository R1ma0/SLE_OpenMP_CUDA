#include <stdio.h>
#include <stdlib.h>

const int BlockSize = 16;

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

__global__ void kernelDown
(
    double *A,
    double *B,
    unsigned int n,
    unsigned int numberRow,
    unsigned int numberColumn
)
{
    int bx = blockIdx.x;
    int by = blockIdx.y;

    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int row = by * BlockSize+ty;
    int column = bx * BlockSize+tx;

    if (numberColumn < column && numberRow < row)
    {
        double glav = A[row * n + numberColumn] / A[numberRow * n + numberColumn];

        if (numberColumn == column)
        {
            B[row] -= B[numberRow] * glav;
        }

        A[row * n + column] -= glav * A[numberRow * n + column];
    }
}

__global__ void kernelUp
(
    double *A,
    double *B,
    unsigned int n,
    unsigned int numberRow,
    unsigned int numberColumn
)
{
    int bx = blockIdx.x;
    int by = blockIdx.y;

    int tx = threadIdx.x;
    int ty = threadIdx.y;

    int row = by * BlockSize+ty;
    int column = bx * BlockSize+tx;

    if (numberColumn == column && numberRow > row)
    {
        double glav = A[row * n + numberColumn] / A[numberRow * n + numberColumn];

        B[row] -= B[numberRow] * glav;
        A[row * n + column] = 0;
    }
}

float solve(double *h_A, double *h_B, double *h_C, unsigned int n)
{
    double *dev_a; 
    double *dev_b;

    cudaMalloc((void **) &dev_a, n * n * sizeof(double));
    cudaMalloc((void **) &dev_b, n * sizeof(double));

    cudaMemcpy(dev_a, h_A, n * n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, h_B, n * sizeof(double), cudaMemcpyHostToDevice);

    dim3 Grid(n / BlockSize, n / BlockSize);
    dim3 Block(BlockSize, BlockSize);

    cudaEvent_t start;
    cudaEvent_t stop;
    float gpuTime = 0;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    for (unsigned int i = 0; i < n - 1; i++)
    {
        kernelDown<<<Grid, Block>>>(dev_a, dev_b, n, i, i);
    }

    for (unsigned int i = n - 1; i >= 0; i--)
    {
        kernelUp<<<Grid,Block>>>(dev_a, dev_b, n, i, i);
    }

    cudaMemcpy(h_A, dev_a, n * n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_B, dev_b, n * sizeof(double), cudaMemcpyDeviceToHost);

    for (unsigned int i = 0; i < n - 1; i++)
    {
        h_C[i] = h_B[i] / h_A[i * n + i];
    }

    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaDeviceReset();

    return gpuTime;
}

int main(int argc, char *argv[])
{
    unsigned int n;

    if (argc <= 1)
    {
        printf("Enter size of the matrix.\n");
        return 1;
    }
    else
    {
        n = atoi(argv[1]);
        printf("Matrix size: %d\n", n);
    }

    // Memory allocation and filling arrays

    double **A = (double **)malloc(n * sizeof(double *));
    for (unsigned int i = 0; i < n; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    double *B = (double *)malloc(n * sizeof(double));
    double *C = (double *)malloc(n * sizeof(double));

    for (unsigned int j = 0; j < n; j++)
    {
        for (unsigned int t = 0; t < n; t++)
        {
            A[j][t] = t;
        }

        B[j] = j;
    }

    // Processing main matrix

    float gpuTime = solve(A, B, C, n);

    printf("GPU time: %.2f milliseconds\n", gpuTime);

    // Freeing memory

    for (unsigned int k = 0; k < n; k++)
    {
        free(A[k]);
    }
    free(A);
    free(B);
    free(C);

    // Checking the correctness of the solution
    
    unsigned int M = 3;
    
    double *correctSolution = (double *)malloc(M * sizeof(double));
    correctSolution[0] = 7.0;
    correctSolution[1] = 5.0;
    correctSolution[2] = 2.0;

    double **a = (double **)malloc(M * sizeof(double *));
    for (unsigned int k = 0; k < M; k++)
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

    // Processing test matrix

    printf("\n=== Checking the correctness of the algorithm  ===\n");

    printf("Test matrix:\n");
    displayMatrix(a, y, M);

    gpuTime = solve(a, y, x, M);

    printf("Resultion solution:\n");
    displaySolution(x, M);

    printf("Correct solution\n");
    displaySolution(correctSolution, M);

    // Freeing memory

    for (unsigned int k = 0; k < M; k++)
    {
        free(a[k]);
    }
    free(a);
    free(y);
    free(x);
    free(correctSolution);

    return 0;
}
