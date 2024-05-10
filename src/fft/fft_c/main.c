#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _WIN64
#include <windows.h>
#endif

#ifdef __linux__
#include <pthread.h>
#endif

#define M_PI 3.14159265358979323846

int fftSize, numLoops;

void fft(int N, double x[N])
{
    double real, imag, tReal, tImag;

    if (N == 2)
    {
        return;
    }

    double *even = malloc(N / 2 * sizeof(double));
    double *odd = malloc(N / 2 * sizeof(double));

    for (int i = 0; i < N / 2; i += 2)
    {
        even[i] = x[2 * i];
        even[i + 1] = x[2 * i + 1];
        odd[i] = x[2 * i + 2];
        odd[i + 1] = x[2 * i + 3];
    }

    fft(N / 2, even);
    fft(N / 2, odd);

    for (int k = 0; k < N / 2; k += 2)
    {
        imag = -2 * M_PI * k / N;
        real = cos(imag);
        imag = sin(imag);

        tReal = real * odd[k] - imag * odd[k + 1];
        tImag = real * odd[k + 1] + imag * odd[k];

        x[k] = even[k] + tReal;
        x[k + 1] = even[k + 1] + tImag;

        x[k + N / 2] = even[k] - tReal;
        x[k + N / 2 + 1] = even[k + 1] - tImag;
    }

    free(even);
    free(odd);
}

void fftLoop(int size, int loops)
{
    srand(time(NULL));

    double *x = malloc(2 * size * sizeof(double));

    for (int loop = 0; loop < loops; loop++)
    {
        for (int i = 0; i < 2 * size; i++)
        {
            x[i] = (double)rand() / (double)RAND_MAX;
        }
        fft(2 * size, x);
    }

    free(x);
}

#ifdef _WIN64
DWORD WINAPI WinThreadProc(LPVOID lpParam)
{
    fftLoop(fftSize, numLoops);
    return 0;
}
#endif

#ifdef __linux__
void *linuxThreadCaller(void *args)
{
    fftLoop(fftSize, numLoops);
    return NULL;
}
#endif

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("Usage: ./fft_c <Input Size> <Loops> <Threads>\n");
        return 1;
    }

    fftSize = 1 << atoi(argv[1]);
    numLoops = atoi(argv[2]);
    int numThreads = atoi(argv[3]);

    if (numThreads == 1)
    {
        fftLoop(fftSize, numLoops);
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
