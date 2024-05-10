class Main {
  static void fft(double[] x) {
    int n = x.length;

    double real, imag, tReal, tImag;

    if (n == 2) {
      return;
    }

    double[] even = new double[n / 2];
    double[] odd = new double[n / 2];

    for (int i = 0; i < n / 2; i += 2) {
      even[i] = x[2 * i];
      even[i + 1] = x[2 * i + 1];
      odd[i] = x[2 * i + 2];
      odd[i + 1] = x[2 * i + 3];
    }

    fft(even);
    fft(odd);

    for (int k = 0; k < n / 2; k += 2) {
      imag = -2 * Math.PI * k / n;
      real = Math.cos(imag);
      imag = Math.sin(imag);

      tReal = real * odd[k] - imag * odd[k + 1];
      tImag = real * odd[k + 1] + imag * odd[k];

      x[k] = even[k] + tReal;
      x[k + 1] = even[k + 1] + tImag;

      x[k + n / 2] = even[k] - tReal;
      x[k + n / 2 + 1] = even[k + 1] - tImag;
    }
  }

  public static void fftLoop(int size, int loops) {
    double[] x = new double[2 * size];

    for (int loop = 0; loop < loops; loop++) {
      for (int i = 0; i < 2 * size; i++) {
        x[i] = Math.random();
      }
      fft(x);
    }
  }

  public static void main(String[] args) {
    if (args.length != 3) {
      System.out.println("Usage: java Main <Input Size> <Loops> <Threads>");
      return;
    }

    int fftSize = 1 << Integer.parseInt(args[0]);
    int numLoops = Integer.parseInt(args[1]);
    int numThreads = Integer.parseInt(args[2]);

    if (numThreads == 1) {
      fftLoop(fftSize, numLoops);
      System.exit(0);
    }

    Thread[] threads = new Thread[numThreads];

    for (int i = 0; i < numThreads; i++) {
      threads[i] = new Thread(() -> fftLoop(fftSize, numLoops));
      threads[i].start();
    }

    for (int i = 0; i < numThreads; i++) {
      try {
        threads[i].join();
      } catch (InterruptedException e) {
        System.out.println("Thread interrupted");
      }
    }
  }
}
