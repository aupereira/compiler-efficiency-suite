import cmath
import math
import multiprocessing
import random
import sys


def rand_complex():
    return complex(random.random(), random.random())


def fft(x):
    N = len(x)

    if N <= 1:
        return x

    even = fft(x[0::2])
    odd = fft(x[1::2])

    T = [cmath.exp(-2j * math.pi * k / N) * odd[k] for k in range(N // 2)]
    return [even[k] + T[k]
            for k in range(N // 2)] + [even[k] - T[k] for k in range(N // 2)]


def fft_loop(size, loops):
    for _ in range(loops):
        x = [rand_complex() for _ in range(size)]
        fft(x)


def main():
    if len(sys.argv) != 4:
        print("Usage: py fft.py <Input Size> <Loops> <Threads>")
        return

    fft_size = 1 << int(sys.argv[1])
    num_loops = int(sys.argv[2])
    num_threads = int(sys.argv[3])

    if num_threads == 1:
        fft_loop(fft_size, num_loops)
        return

    processes = []

    for i in range(num_threads):
        processes.append(multiprocessing.Process(
            target=fft_loop, args=(fft_size, num_loops)
        ))
        processes[i].start()

    for process in processes:
        process.join()


if __name__ == '__main__':
    main()
