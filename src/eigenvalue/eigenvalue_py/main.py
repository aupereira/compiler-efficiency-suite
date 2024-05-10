import math
import multiprocessing
import random
import sys


def random_symmetric(n):
    A = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        A[i][i] = random.random()
        for j in range(i + 1, n):
            A[i][j] = random.random()
            A[j][i] = A[i][j]

    return A


def sgn(x):
    if x >= 0:
        return 1
    else:
        return -1


def norm2(v):
    norm = 0
    for i in v:
        norm += i * i
    return math.sqrt(norm)


def v_div_s(v, s):
    return [i / s for i in v]


def S2OPAI(v):
    length = len(v)

    A = [[-2 * v[i] * v[j] for i in range(length)] for j in range(length)]

    for i in range(length):
        A[i][i] += 1

    return A


def transpose(A):
    return [list(i) for i in zip(*A)]


def matrix_mul(A, B):
    BTMP = transpose(B)
    return [[sum(a * b for a, b in zip(row_a, col_b))
             for col_b in BTMP] for row_a in A]


def tridiagonalize(A):
    N = len(A)

    v = [0 for _ in range(N)]

    for k in range(N - 2):
        v[0:k + 1] = [0 for _ in range(k + 1)]
        v[k + 1:N] = [A[j][k] for j in range(k + 1, N)]

        alpha = -sgn(v[k + 1]) * norm2(v)
        r = 2 * math.sqrt(0.5 * ((alpha * alpha) - (alpha * v[k + 1])))

        v[k + 1] = (v[k + 1] - alpha) / r
        v[k + 2:N] = [i / r for i in v[k + 2:N]]

        P = S2OPAI(v)
        A = matrix_mul(matrix_mul(P, A), P)

    return A


def householder(A):
    N = len(A)

    v = [0 for _ in range(N)]

    s = A[N - 1][N - 1]
    t = A[N - 2][N - 2]
    x = A[N - 2][N - 1]
    d = (t - s) / 2
    shift = s - x * x / (d + sgn(d) * math.sqrt(d * d + x * x))

    for i in range(N):
        A[i][i] -= shift

    for k in range(N - 1):
        v[0:k] = [0 for _ in range(k)]
        v[k:N] = [A[j][k] for j in range(k, N)]

        v[k] += norm2(v) * sgn(v[k + 1])
        norm = norm2(v)
        if norm != 0:
            v = v_div_s(v, norm2(v))
        H = S2OPAI(v)
        A = matrix_mul(H, A)
        if (k == 0):
            Q = H
        else:
            Q = matrix_mul(Q, H)

    A = matrix_mul(A, Q)

    for i in range(N):
        A[i][i] += shift

    return A


def eig_QR(iter, A):
    vals = [0 for _ in range(len(A))]

    for ind in range(len(A), 1, -1):
        for _ in range(iter):
            A = householder(A)
        vals[ind - 1] = A[ind - 1][ind - 1]
        A = [row[:-1] for row in A[:-1]]

    vals[0] = A[0][0]

    return vals


def eigen_loop(size, iter, loops):
    for _ in range(loops):
        eig_QR(iter, tridiagonalize(random_symmetric(size)))


def main():
    if len(sys.argv) != 5:
        print("Usage: py main.py <Matrix Size> <Convergence Iterations> <Loops> <Threads>")
        return

    m_size = int(sys.argv[1])
    c_iter = int(sys.argv[2])
    num_loops = int(sys.argv[3])
    num_threads = int(sys.argv[4])

    if num_threads == 1:
        eigen_loop(m_size, c_iter, num_loops)
        return

    processes = []

    for i in range(num_threads):
        processes.append(multiprocessing.Process(
            target=eigen_loop, args=(m_size, c_iter, num_loops)
        ))
        processes[i].start()

    for process in processes:
        process.join()


if __name__ == '__main__':
    main()
