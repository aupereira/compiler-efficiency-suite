using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

class Program
{

    private static Random rnd = new();

    static void MCopy(int len, double[,] FROM, double[,] TO)
    {
        for (int i = 0; i < len; i++)
        {
            for (int j = 0; j < len; j++)
            {
                TO[i, j] = FROM[i, j];
            }
        }
    }

    static double Sgn(double x) => x >= 0.0 ? 1.0 : -1.0;

    static void RandomSymmetric(double[,] A)
    {
        int len = A.GetLength(0);

        for (int i = 0; i < len; i++)
        {
            A[i, i] = rnd.NextDouble();
            for (int j = i + 1; j < len; j++)
            {
                A[i, j] = rnd.NextDouble();
                A[j, i] = A[i, j];
            }
        }
    }

    static void MatrixMultiply(
        int len,
        double[,] A,
        double[,] B,
        double[,] C,
        double[,] BTMP
    ) {
        double sum;

        for (int i = 0; i < len; i++)
        {
            for (int j = 0; j < len; j++)
            {
                BTMP[i, j] = B[j, i];
            }
        }

        for (int i = 0; i < len; i++)
        {
            for (int j = 0; j < len; j++)
            {
                sum = 0;
                for (int k = 0; k < len; k++)
                {
                    sum += A[i, k] * BTMP[j, k];
                }
                C[i, j] = sum;
            }
        }
    }

    static double Norm2(int len, double[] v)
    {
        double norm = 0;
        for (int i = 0; i < len; i++)
        {
            norm += v[i] * v[i];
        }
        return Math.Sqrt(norm);
    }

    static void S2OPAI(int len, double[] v, double[,] A)
    {
        for (int i = 0; i < len; i++)
        {
            for (int j = 0; j < len; j++)
            {
                A[i, j] = -2.0 * v[i] * v[j];
            }
        }

        for (int i = 0; i < len; i++)
        {
            A[i, i] += 1.0;
        }
    }

    static void VDivS(int len, double[] v, double s)
    {
        for (int i = 0; i < len; i++)
        {
            v[i] /= s;
        }
    }

    static void Tridiagonalize(
        int N,
        double[,] A,
        double[,] P,
        double[,] TEMP,
        double[,] TPOSE,
        double[] v
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
                v[j] = A[j, k];
            }

            alpha = -Sgn(v[k + 1]) * Norm2(N, v);
            r = 2.0 * Math.Sqrt(0.5 * ((alpha * alpha) - (alpha * v[k + 1])));

            v[k + 1] = (v[k + 1] - alpha) / r;
            for (int j = k + 2; j < N; j++)
            {
                v[j] /= r;
            }

            S2OPAI(N, v, P);
            MatrixMultiply(N, P, A, TEMP, TPOSE);
            MatrixMultiply(N, TEMP, P, A, TPOSE);
        }
    }

    static void Householder(
        int len,
        double[,] A,
        double[,] A1,
        double[,] Q,
        double[,] H,
        double[,] TEMP,
        double[,] TPOSE,
        double[] v
    ) {
        MCopy(len, A, A1);
        double[,] SWAP;

        double s = A1[len - 1, len - 1];
        double t = A1[len - 2, len - 2];
        double x = A1[len - 2, len - 1];
        double d = (t - s) / 2.0;
        double shift = s - x * x / (d + Sgn(d) * Math.Sqrt(d * d + x * x));

        for (int i = 0; i < len; i++)
        {
            A1[i, i] -= shift;
        }

        for (int k = 0; k < len - 1; k++)
        {
            for (int j = 0; j < k; j++)
            {
                v[j] = 0.0;
            }
            for (int j = k; j < len; j++)
            {
                v[j] = A1[j, k];
            }

            v[k] += Norm2(len, v) * Sgn(v[k + 1]);
            VDivS(len, v, Norm2(len, v));

            if (k == 0)
            {
                S2OPAI(len, v, Q);
                MatrixMultiply(len, Q, A1, TEMP, TPOSE);
                SWAP = A1;
                A1 = TEMP;
                TEMP = SWAP;
            }
            else
            {
                S2OPAI(len, v, H);
                MatrixMultiply(len, H, A1, TEMP, TPOSE);
                SWAP = A1;
                A1 = TEMP;
                TEMP = SWAP;
                MatrixMultiply(len, Q, H, TEMP, TPOSE);
                SWAP = Q;
                Q = TEMP;
                TEMP = SWAP;
            }
        }

        MatrixMultiply(len, A1, Q, A, TPOSE);

        for (int i = 0; i < len; i++)
        {
            A[i, i] += shift;
        }
    }

    static void EigQR(
        int iter,
        double[,] A,
        double[] eigVals,
        double[,] T1,
        double[,] T2,
        double[,] T3,
        double[,] T4,
        double[,] T5,
        double[] v
    ) {
        int N = A.GetLength(0);

        for (int ind = N; ind > 1; ind--)
        {
            for (int i = 0; i < iter; i++)
            {
                Householder(ind, A, T1, T2, T3, T4, T5, v);
            }
            eigVals[ind - 1] = A[ind - 1, ind - 1];
        }

        eigVals[0] = A[0, 0];
    }

    static void EigenLoop(int size, int iter, int loops)
    {
        double[,] A = new double[size, size];

        double[,] T1 = new double[size, size];
        double[,] T2 = new double[size, size];
        double[,] T3 = new double[size, size];
        double[,] T4 = new double[size, size];
        double[,] T5 = new double[size, size];
        double[] v = new double[size];

        double[] eigvals = new double[size];

        for (int i = 0; i < loops; i++)
        {
            RandomSymmetric(A);
            Tridiagonalize(size, A, T1, T2, T3, v);
            EigQR(iter, A, eigvals, T1, T2, T3, T4, T5, v);
        }
    }

    static void Main(string[] args)
    {
        if (args.Length != 4)
        {
            Console.WriteLine("Usage: ./eigenvalue_cs <Matrix Size> <Convergence Iterations> <Loops> <Threads>");
            return;
        }

        int mSize = int.Parse(args[0]);
        int cIter = int.Parse(args[1]);
        int numLoops = int.Parse(args[2]);
        int numThreads = int.Parse(args[3]);

        if (numThreads == 1)
        {
            EigenLoop(mSize, cIter, numLoops);
            Environment.Exit(0);
        }

        Thread[] threads = new Thread[numThreads];

        for (int i = 0; i < numThreads; i++)
        {
            threads[i] = new Thread(() => EigenLoop(mSize, cIter, numLoops));
            threads[i].Start();
        }

        for (int i = 0; i < numThreads; i++)
        {
            threads[i].Join();
        }
    }
}
