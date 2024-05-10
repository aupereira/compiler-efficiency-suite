#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#define sgn(x) ((x) >= 0.0 ? 1.0 : -1.0)

typedef std::vector<double> VectorDouble;
typedef std::vector<std::vector<double>> MatrixDouble;

void rndSymmetric(MatrixDouble &A)
{
    for (auto i = 0; i < A.size(); i++)
    {
        A[i][i] = rand() / (double)RAND_MAX;
        for (auto j = i + 1; j < A.size(); j++)
        {
            A[i][j] = rand() / (double)RAND_MAX;
            A[j][i] = A[i][j];
        }
    }
}

double norm2(auto len, VectorDouble &v)
{
    double norm = 0;
    for (auto i = 0; i < len; i++)
    {
        norm += v[i] * v[i];
    }
    return std::sqrt(norm);
}

void vDivS(auto len, VectorDouble &v, double s)
{
    for (auto i = 0; i < len; i++)
    {
        v[i] /= s;
    }
}

void S2OPAI(auto len, VectorDouble &v, MatrixDouble &A)
{
    for (auto i = 0; i < len; i++)
    {
        for (auto j = 0; j < len; j++)
        {
            A[i][j] = -2.0 * v[i] * v[j];
        }
    }

    for (auto i = 0; i < len; i++)
    {
        A[i][i] += 1.0;
    }
}

void mMul(
    auto len,
    MatrixDouble &A,
    MatrixDouble &B,
    MatrixDouble &C,
    MatrixDouble &BTMP
) {
    double sum;

    for (auto i = 0; i < len; i++)
    {
        for (auto j = 0; j < len; j++)
        {
            BTMP[i][j] = B[j][i];
        }
    }

    for (auto i = 0; i < len; i++)
    {
        for (auto j = 0; j < len; j++)
        {
            sum = 0;
            for (auto k = 0; k < len; k++)
            {
                sum += A[i][k] * BTMP[j][k];
            }
            C[i][j] = sum;
        }
    }
}

void tridiagonalize(
    MatrixDouble &A,
    MatrixDouble &P,
    MatrixDouble &TEMP,
    MatrixDouble &TPOSE,
    VectorDouble &v
) {
    auto N = A.size();

    double alpha;
    double r;

    for (auto k = 0; k < N - 2; k++)
    {
        for (auto j = 0; j < k + 1; j++)
        {
            v[j] = 0;
        }
        for (auto j = k + 1; j < N; j++)
        {
            v[j] = A[j][k];
        }

        alpha = -sgn(v[k + 1]) * norm2(N, v);
        r = 2 * std::sqrt(0.5 * ((alpha * alpha) - (alpha * v[k + 1])));

        v[k + 1] = (v[k + 1] - alpha) / r;
        for (auto j = k + 2; j < N; j++)
        {
            v[j] /= r;
        }

        S2OPAI(N, v, P);
        mMul(N, P, A, TEMP, TPOSE);
        mMul(N, TEMP, P, A, TPOSE);
    }
}

void householder(
    auto len,
    MatrixDouble &A,
    MatrixDouble *A1,
    MatrixDouble *Q,
    MatrixDouble *H,
    MatrixDouble *TEMP,
    MatrixDouble &TPOSE,
    VectorDouble &v
) {
    *A1 = A;
    MatrixDouble *SWAP;

    double s = (*A1)[len - 1][len - 1];
    double t = (*A1)[len - 2][len - 2];
    double x = (*A1)[len - 2][len - 1];
    double d = (t - s) / 2;
    double shift = s - x * x / (d + sgn(d) * std::sqrt(d * d + x * x));

    for (auto i = 0; i < len; i++)
    {
        (*A1)[i][i] -= shift;
    }

    for (auto k = 0; k < len - 1; k++)
    {
        for (auto j = 0; j < k; j++)
        {
            v[j] = 0;
        }
        for (auto j = k; j < len; j++)
        {
            v[j] = (*A1)[j][k];
        }

        v[k] += norm2(len, v) * sgn(v[k + 1]);
        vDivS(len, v, norm2(len, v));

        if (k == 0)
        {
            S2OPAI(len, v, *Q);
            mMul(len, *Q, *A1, *TEMP, TPOSE);
            SWAP = A1;
            A1 = TEMP;
            TEMP = SWAP;
        }
        else
        {
            S2OPAI(len, v, *H);
            mMul(len, *H, *A1, *TEMP, TPOSE);
            SWAP = A1;
            A1 = TEMP;
            TEMP = SWAP;
            mMul(len, *Q, *H, *TEMP, TPOSE);
            SWAP = Q;
            Q = TEMP;
            TEMP = SWAP;
        }
    }

    mMul(len, *A1, *Q, A, TPOSE);

    for (auto i = 0; i < len; i++)
    {
        A[i][i] += shift;
    }
}

void eigQR(
    int iter,
    MatrixDouble &A,
    VectorDouble &eigVals,
    MatrixDouble &T1,
    MatrixDouble &T2,
    MatrixDouble &T3,
    MatrixDouble &T4,
    MatrixDouble &T5,
    VectorDouble &v
) {
    for (int ind = A.size(); ind > 1; ind--)
    {
        for (int i = 0; i < iter; i++)
        {
            householder(ind, A, &T1, &T2, &T3, &T4, T5, v);
        }
        eigVals[ind - 1] = A[ind - 1][ind - 1];
    }

    eigVals[0] = A[0][0];
}

void eigenLoop(int size, int iter, int loops)
{
    srand(time(0));

    MatrixDouble A(size, VectorDouble(size));

    MatrixDouble T1(size, VectorDouble(size));
    MatrixDouble T2(size, VectorDouble(size));
    MatrixDouble T3(size, VectorDouble(size));
    MatrixDouble T4(size, VectorDouble(size));
    MatrixDouble T5(size, VectorDouble(size));
    VectorDouble v(size);

    VectorDouble eigvals(size);

    for (int loop = 0; loop < loops; loop++)
    {
        rndSymmetric(A);
        tridiagonalize(A, T1, T2, T3, v);
        eigQR(iter, A, eigvals, T1, T2, T3, T4, T5, v);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cout << "Usage: ./eigenvalue_cpp <Matrix Size> <Convergence Iterations> <Loops> <Threads>" << std::endl;
        exit(1);
    }

    int mSize = std::stoi(argv[1]);
    int cIter = std::stoi(argv[2]);
    int numLoops = std::stoi(argv[3]);
    int numThreads = std::stoi(argv[4]);

    if (numThreads == 1)
    {
        eigenLoop(mSize, cIter, numLoops);
        exit(0);
    }

    std::vector<std::thread> threads(numThreads);

    for (auto i = 0; i < numThreads; i++)
    {
        threads[i] = std::thread([mSize, cIter, numLoops]()
                                 { eigenLoop(mSize, cIter, numLoops); });
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    exit(0);
}
