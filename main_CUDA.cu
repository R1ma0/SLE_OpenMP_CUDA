#include "utils.h"

__device__ void findMaxRowValue(
    double *max, double *A, unsigned int *N, unsigned int *index, 
    unsigned int *k
)
{
    unsigned int xIdx = threadIdx.x + k + 1;

    if (xIdx < N)
    {
        if (abs(A[xIdx][k]) > max)
        {
            max = abs(A[xIdx][k]);
            index = xIdx;
        }
    }
}

__device__ void rearrangingStrings(
    double *A, unsigned int *k, unsigned int *N, unsigned int *index
)
{
    double tmp;
    unsigned int xIdx = threadIdx.x;

    if (xIdx < N)
    {
        tmp = A[k][xIdx];
        A[k][xIdx] = A[index][xIdx];
        A[index][xIdx] = tmp;
    }
}

__device__ void normalize(
    double *A, double *Y, double *tmp, unsigned int *k, unsigned int *N
)
{
    unsigned int xIdx = threadIdx.x + k;
    unsigned int yIdx = threadIdx.y + k;

    if (xIdx < N)
    {
        tmp = A[xIdx][k];
        
        if (abs(tmp) < DBL_EPSILON)
        {
            continue;
        }

        if (yIdx < N)
        {
            A[xIdx][yIdx] = A[xIdx][yIdx] / tmp;
        }
        
        Y[xIdx] = Y[xIdx] / tmp;

        if (xIdx == k)
        {
            continue;
        }

        if (yIdx < N)
        {
            A[xIdx][yIdx] = A[xIdx][yIdx] - A[k][yIdx];
        }

        Y[xIdx] = Y[xIdx] - Y[k];
    }
}

__global__ void useGaussMethod(
    double *A, double *Y, double *X, unsigned int *N, double *max,
    unsigned int *k, unsigned int *index, dim3 *threads, dim3 *blocks,
    double *tmp
)
{
    unsigned int i;
    unsigned int j;

    while (k < N)
    {
        // Search max A[i][k]

        max = abs(A[k][k]);
        index = k;

        findMaxRowValue<<<blocks, threads>>>(max, A, N, index, k);

        // Rearranging strings

        rearrangingStrings<<<blocks, threads>>>(A, k, N, index);

        tmp = Y[k];
        Y[k] = Y[index];
        Y[index] = tmp;

        // Normalization

        normalize<<<blocks, threads>>>(A, Y, tmp, k, N);

        k++;
    }

    // Reverse substitution

    unsigned int t;
    unsigned int g;

    for (t = N - 1; t >= 0; t--)
    {
        X[t] = Y[t];
        
        for (g = 0; g < t; g++)
        {
            Y[g] -= A[g][t] * X[t];
        }
    }
}

void useGaussMethod(double **A, double *Y, double *X, unsigned int N)
{
    double max,
           tmp;
    int g,
        t;
    unsigned int k = 0,
                 i,
                 j,
                 index;

    while (k < N)
    {
        // Search max A[i][k]

        max = abs(A[k][k]);
        index = k;

#pragma omp parallel for schedule(dynamic, 32) private(i)
        for (i = k + 1; i < N; i++)
        {
            if (abs(A[i][k]) > max)
            {
                max = abs(A[i][k]);
                index = i;
            }
        }

        // Rearranging strings

#pragma omp parallel for schedule(dynamic, 32) private(j, tmp)
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

#pragma omp parallel for schedule(dynamic, 32) private(i, j, tmp)
        for (i = k; i < N; i++)
        {
            tmp = A[i][k];
            
            if (abs(tmp) < DBL_EPSILON)
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

    for (t = N - 1; t >= 0; t--)
    {
        X[t] = Y[t];

#pragma omp parallel for schedule(dynamic, 32) private(g)
        for (g = 0; g < t; g++)
        {
            Y[g] -= A[g][t] * X[t];
        }
    }
}

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

    // CUDA

    double gpuTime = 0.0;
    cudaEvent_t start;
    cudaEvent_t stop;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    dim3 threads = dim3(512, 1);
    dim3 blocks = dim3(N / threads.x, 1);

    dim3 *threadsDev;
    dim3 *blocksDev;

    double *aDev = NULL;
    double *xDev = NULL;
    double *yDev = NULL;

    double max = 0.0;
    double tmp = 0.0;
    unsigned int k = 0;
    unsigned int index = 0;

    double *maxDev;
    double *tmpDev;
    unsigned int *kDev;
    unsigned int *indexDev;

    cudaMalloc((void **)&threadsDev, sizeof(dim3));
    cudaMalloc((void **)&blocksDev, sizeof(dim3));
    cudaMalloc((void **)&aDev, N * N * sizeof(double));
    cudaMalloc((void **)&xDev, N * sizeof(double));
    cudaMalloc((void **)&yDev, N * sizeof(double));
    cudaMalloc((void **)&NDev, sizeof(unsigned int));
    cudaMalloc((void **)&maxDev, sizeof(double));
    cudaMalloc((void **)&tmpDev, sizeof(double));
    cudaMalloc((void **)&kDev, sizeof(unsigned int));
    cudaMalloc((void **)&indexDev, sizeof(unsigned int));

    cudaMemcpy(threadsDev, threads, sizeof(dim3), cudaMemcpyHostToDevice);
    cudaMemcpy(blocksDev, blocks, sizeof(dim3), cudaMemcpyHostToDevice);
    cudaMemcpy(aDev, A, N * N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(xDev, X, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(yDev, Y, N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(NDev, N, sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(maxDev, max, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(tmpDev, tmp, sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(kDev, k, sizeof(unsigned int), cudaMemcpyHostToDevice);
    cudaMemcpy(indexDev, sizeof(unsigned int), cudaMemcpyHostToDevice);

    cudaEventRecord(start, 0);
    useGaussMethod<<<blocks, threads>>>(
        aDev, yDev, xDev, NDev, maxDev, kDev, indexDev, threadsDev, blocksDev,
        tmpDev
    );
    cudaEventRecord(stop, 0);
    
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    printf("GPU time: %.2f milliseconds\n", gpuTime);

    cudaMemcpy(Y, yDev, N * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(aDev);
    cudaFree(xDev);
    cudaFree(yDev);
    cudaFree(NDev);
    cudaFree(maxDev);
    cudaFree(tmpDev);
    cudaFree(kDev);
    cudaFree(indexDev);
    cudaFree(threadsDev);
    cudaFree(blocksDev);

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

    for (k = 0; k < M; k++)
    {
        free(a[k]);
    }
    free(a);
    free(y);
    free(x);
    free(correctSolution);

    for (i = 0; i < N; i++)
    {
        free(A[i]);
    }
    free(A);
    free(Y);
    free(X);

    return 0;
}

