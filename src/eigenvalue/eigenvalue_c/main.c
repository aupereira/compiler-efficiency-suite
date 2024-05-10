#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN64
#include <Windows.h>
#endif

#ifdef __linux__
#include <pthread.h>
#endif

#define sgn(x) (x >= 0.0 ? 1.0 : -1.0)

int mSize, cIter, numLoops;

void mCopy(int N, int len, double FROM[N][N], double TO[N][N])
{
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            TO[i][j] = FROM[i][j];
        }
    }
}

void rndSymmetric(int N, double A[N][N])
{
    for (int i = 0; i < N; i++)
    {
        A[i][i] = rand() / (double)RAND_MAX;
        for (int j = i + 1; j < N; j++)
        {
            A[i][j] = rand() / (double)RAND_MAX;
            A[j][i] = A[i][j];
        }
    }
}

double norm2(int N, int len, double v[N])
{
    double norm = 0.0;
    for (int i = 0; i < len; i++)
    {
        norm += v[i] * v[i];
    }
    return sqrt(norm);
}

void vDivS(int N, int len, double v[N], double s)
{
    for (int i = 0; i < len; i++)
    {
        v[i] /= s;
    }
}

void S2OPAI(int N, int len, double v[N], double A[N][N])
{
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            A[i][j] = -2.0 * v[i] * v[j];
        }
    }

    for (int i = 0; i < len; i++)
    {
        A[i][i] += 1.0;
    }
}

void mMul(
    int N,
    int len,
    double A[N][N],
    double B[N][N],
    double C[N][N],
    double BTMP[N][N]
) {
    double sum;

    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            BTMP[i][j] = B[j][i];
        }
    }

    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            sum = 0.0;
            for (int k = 0; k < len; k++)
            {
                sum += A[i][k] * BTMP[j][k];
            }
            C[i][j] = sum;
        }
    }
}

void tridiagonalize(
    int N,
    double A[N][N],
    double P[N][N],
    double TEMP[N][N],
    double TPOSE[N][N],
    double v[N]
) {
    double alpha;
    double r;

    for (int k = 0; k < N - 2; k++)
    {
        for (int j = 0; j < k + 1; j++)
        {
            v[j] = 0.0;
        }
        for (int j = k + 1; j < N; j++)
        {
            v[j] = A[j][k];
        }

        alpha = -sgn(v[k + 1]) * norm2(N, N, v);
        r = 2.0 * sqrt(0.5 * ((alpha * alpha) - (alpha * v[k + 1])));

        v[k + 1] = (v[k + 1] - alpha) / r;
        for (int j = k + 2; j < N; j++)
        {
            v[j] /= r;
        }

        S2OPAI(N, N, v, P);
        mMul(N, N, P, A, TEMP, TPOSE);
        mMul(N, N, TEMP, P, A, TPOSE);
    }
}

void householder(
    int N,
    int len,
    double A[N][N],
    double A1[N][N],
    double Q[N][N],
    double H[N][N],
    double TEMP[N][N],
    double TPOSE[N][N],
    double v[N]
) {
    mCopy(N, len, A, A1);
    double(*SWAP)[N];

    double s = A1[len - 1][len - 1];
    double t = A1[len - 2][len - 2];
    double x = A1[len - 2][len - 1];
    double d = (t - s) / 2.0;
    double shift = s - x * x / (d + sgn(d) * sqrt(d * d + x * x));

    for (int i = 0; i < len; i++)
    {
        A1[i][i] -= shift;
    }

    for (int k = 0; k < len - 1; k++)
    {
        for (int j = 0; j < k; j++)
        {
            v[j] = 0.0;
        }
        for (int j = k; j < len; j++)
        {
            v[j] = A1[j][k];
        }

        v[k] += norm2(N, len, v) * sgn(v[k + 1]);
        vDivS(N, len, v, norm2(N, len, v));

        if (k == 0)
        {
            S2OPAI(N, len, v, Q);
            mMul(N, len, Q, A1, TEMP, TPOSE);
            SWAP = A1;
            A1 = TEMP;
            TEMP = SWAP;
        }
        else
        {
            S2OPAI(N, len, v, H);
            mMul(N, len, H, A1, TEMP, TPOSE);
            SWAP = A1;
            A1 = TEMP;
            TEMP = SWAP;
            mMul(N, len, Q, H, TEMP, TPOSE);
            SWAP = Q;
            Q = TEMP;
            TEMP = SWAP;
        }
    }

    mMul(N, len, A1, Q, A, TPOSE);

    for (int i = 0; i < len; i++)
    {
        A[i][i] += shift;
    }
}

void eigQR(
    int N,
    int iter,
    double A[N][N],
    double eigvals[N],
    double T1[N][N],
    double T2[N][N],
    double T3[N][N],
    double T4[N][N],
    double T5[N][N],
    double v[N]
) {
    for (int ind = N; ind > 1; ind--)
    {
        for (int i = 0; i < iter; i++)
        {
            householder(N, ind, A, T1, T2, T3, T4, T5, v);
        }
        eigvals[ind - 1] = A[ind - 1][ind - 1];
    }

    eigvals[0] = A[0][0];
}

void eigenLoop(int size, int iter, int numLoops)
{
    srand(time(NULL));

    double(*A)[size] = malloc(mSize * sizeof *A);

    double(*T1)[size] = malloc(size * sizeof *T1);
    double(*T2)[size] = malloc(size * sizeof *T2);
    double(*T3)[size] = malloc(size * sizeof *T3);
    double(*T4)[size] = malloc(size * sizeof *T4);
    double(*T5)[size] = malloc(size * sizeof *T5);
    double *v = (double *)malloc(size * sizeof(double));

    double *eigvals = (double *)malloc(size * sizeof(double));

    for (size_t i = 0; i < numLoops; i++)
    {
        rndSymmetric(size, A);
        tridiagonalize(size, A, T1, T2, T3, v);
        eigQR(size, iter, A, eigvals, T1, T2, T3, T4, T5, v);
    }

    free(T1);
    free(T2);
    free(T3);
    free(T4);
    free(v);
    free(eigvals);
}

#ifdef _WIN64
DWORD WINAPI WinThreadProc()
{
    eigenLoop(mSize, cIter, numLoops);
    return 0;
}
#endif

#ifdef __linux__
void *linuxThreadCaller(void *args)
{
    eigenLoop(mSize, cIter, numLoops);
    return NULL;
}
#endif

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("Usage: ./eigenvalue_c <Matrix Size> <Convergence Iterations> <Loops> <Threads>\n");
        return 1;
    }

    mSize = atoi(argv[1]);
    cIter = atoi(argv[2]);
    numLoops = atoi(argv[3]);
    int numThreads = atoi(argv[4]);

    if (numThreads == 1)
    {
        eigenLoop(mSize, cIter, numLoops);
        return 0;
    }

#ifdef _WIN64
    HANDLE *hThreads = malloc(numThreads * sizeof(HANDLE));

    for (int i = 0; i < numThreads; i++)
    {
        hThreads[i] = CreateThread(
            NULL,
            0,
            WinThreadProc,
            NULL,
            0,
            NULL);

        if (hThreads[i] == NULL)
        {
            printf("CreateThread failed (%lu)\n", GetLastError());
            return 1;
        }
    }

    WaitForMultipleObjects(numThreads, hThreads, TRUE, INFINITE);

    for (int i = 0; i < numThreads; i++)
    {
        CloseHandle(hThreads[i]);
    }

    free(hThreads);
#endif

#ifdef __linux__
    pthread_t *threads = malloc(numThreads * sizeof(pthread_t));

    for (int i = 0; i < numThreads; i++)
    {
        pthread_create(&threads[i], NULL, linuxThreadCaller, NULL);
    }

    for (int i = 0; i < numThreads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    free(threads);
#endif

    return 0;
}
