using System;
using System.Threading;

class Program
{
    static void FFT(ref double[] x)
    {
        double real, imag, tReal, tImag;

        int n = x.Length;

        if (n == 2)
        {
            return;
        }

        double[] even = new double[n / 2];
        double[] odd = new double[n / 2];

        for (int i = 0; i < n / 2; i += 2)
        {
            even[i] = x[2 * i];
            even[i + 1] = x[2 * i + 1];
            odd[i] = x[2 * i + 2];
            odd[i + 1] = x[2 * i + 3];
        }

        FFT(ref even);
        FFT(ref odd);

        for (int k = 0; k < n / 2; k += 2)
        {
            imag = -2 * Math.PI * k / n;
            real = Math.Cos(imag);
            imag = Math.Sin(imag);

            tReal = real * odd[k] - imag * odd[k + 1];
            tImag = real * odd[k + 1] + imag * odd[k];

            x[k] = even[k] + tReal;
            x[k + 1] = even[k + 1] + tImag;

            x[k + n / 2] = even[k] - tReal;
            x[k + n / 2 + 1] = even[k + 1] - tImag;
        }
    }

    static void FFTLoop(int size, int loops)
    {
        double[] x = new double[2 * size];
        Random rand = new Random();

        for (int loop = 0; loop < loops; loop++)
        {
            for (int i = 0; i < 2 * size; i++)
            {
                x[i] = rand.NextDouble();
            }
            FFT(ref x);
        }
        return;
    }

    static void Main(string[] args)
    {
        if (args.Length != 3)
        {
            Console.WriteLine("Usage: ./fft_cs <Input Size> <Loops> <Threads>");
            return;
        }

        int fftSize = 1 << int.Parse(args[0]);
        int numLoops = int.Parse(args[1]);
        int numThreads = int.Parse(args[2]);

        if (numThreads == 1)
        {
            FFTLoop(fftSize, numLoops);
            Environment.Exit(0);
        }

        Thread[] threads = new Thread[numThreads];

        for (int i = 0; i < numThreads; i++)
        {
            threads[i] = new Thread(() => FFTLoop(fftSize, numLoops));
            threads[i].Start();
        }

        foreach (Thread thread in threads)
        {
            thread.Join();
        }
    }
}
