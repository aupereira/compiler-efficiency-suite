#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#define M_PI 3.14159265358979323846

void fft(std::vector<double> &x)
{
    int n = x.size();

    double real, imag, tReal, tImag;

    if (n == 2)
    {
        return;
    }

    std::vector<double> even(n / 2);
    std::vector<double> odd(n / 2);

    for (int i = 0; i < n / 2; i += 2)
    {
        even[i] = x[2 * i];
        even[i + 1] = x[2 * i + 1];
        odd[i] = x[2 * i + 2];
        odd[i + 1] = x[2 * i + 3];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < n / 2; k += 2)
    {
        imag = -2 * M_PI * k / n;
        real = cos(imag);
        imag = sin(imag);

        tReal = real * odd[k] - imag * odd[k + 1];
        tImag = real * odd[k + 1] + imag * odd[k];

        x[k] = even[k] + tReal;
        x[k + 1] = even[k + 1] + tImag;

        x[k + n / 2] = even[k] - tReal;
        x[k + n / 2 + 1] = even[k + 1] - tImag;
    }
}

void fftLoop(int size, int loops)
{
    srand(time(0));

    std::vector<double> x(2 * size);

    for (int loop = 0; loop < loops; loop++)
    {
        for (int i = 0; i < 2 * size; i++)
        {
            x[i] = (double)rand() / (double)RAND_MAX;
        }
        fft(x);
    }
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cout << "Usage: ./fft_cpp <Input Size> <Loops> <Threads>" << std::endl;
        exit(1);
    }

    auto fftSize = 1 << std::stoi(argv[1]);
    auto numLoops = std::stoi(argv[2]);
    auto numThreads = std::stoi(argv[3]);

    if (numThreads == 1)
    {
        fftLoop(fftSize, numLoops);
        exit(0);
    }

    std::vector<std::thread> threads(numThreads);

    for (std::thread &thread : threads)
    {
        thread = std::thread(fftLoop, fftSize, numLoops);
    }

    for (std::thread &thread : threads)
    {
        thread.join();
    }

    exit(0);
}
