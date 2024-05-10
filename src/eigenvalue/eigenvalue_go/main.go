package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"sync"
)

func create2DSlice(size int) [][]float64 {
	matrix := make([][]float64, size)
	rows := make([]float64, size*size)
	for i := 0; i < size; i++ {
		matrix[i] = rows[i*size : (i+1)*size]
	}
	return matrix
}

func copy2DSlice(dest *[][]float64, src *[][]float64) {
	for i := 0; i < len(*src); i++ {
		copy((*dest)[i], (*src)[i])
	}
}

func sgn(x float64) float64 {
	if x >= 0 {
		return 1
	} else {
		return -1
	}
}

func randSymmetric(A *[][]float64) {
	n := len(*A)

	for i := 0; i < n; i++ {
		(*A)[i][i] = rand.Float64()
		for j := 0; j < i; j++ {
			(*A)[i][j] = rand.Float64()
			(*A)[j][i] = (*A)[i][j]
		}
	}
}

func norm2(len int, v *[]float64) float64 {
	sum := 0.0

	for i := 0; i < len; i++ {
		sum += (*v)[i] * (*v)[i]
	}

	return math.Sqrt(sum)
}

func S2OPAI(len int, v *[]float64, A *[][]float64) {
	for i := 0; i < len; i++ {
		for j := 0; j < len; j++ {
			(*A)[i][j] = -2 * (*v)[i] * (*v)[j]
		}
	}

	for i := 0; i < len; i++ {
		(*A)[i][i] += 1
	}
}

func vDivS(len int, v *[]float64, s float64) {
	for i := 0; i < len; i++ {
		(*v)[i] /= s
	}
}

func matMul(
	n int,
	A *[][]float64,
	B *[][]float64,
	C *[][]float64,
	BTMP *[][]float64,
) {
	var sum float64

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			(*BTMP)[i][j] = (*B)[j][i]
		}
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			sum = 0.0
			for k := 0; k < n; k++ {
				sum += (*A)[i][k] * (*BTMP)[j][k]
			}
			(*C)[i][j] = sum
		}
	}
}

func tridiagonalize(
	n int,
	A *[][]float64,
	P *[][]float64,
	TEMP *[][]float64,
	TPOSE *[][]float64,
	v *[]float64,
) {
	var alpha float64
	var r float64

	for k := 0; k < n-2; k++ {
		for j := 0; j < n; j++ {
			(*v)[j] = 0
		}
		for j := k + 1; j < n; j++ {
			(*v)[j] = (*A)[j][k]
		}

		alpha = -sgn((*v)[k+1]) * norm2(n, v)
		r = 2 * math.Sqrt(0.5*((alpha*alpha)-(alpha*(*v)[k+1])))

		(*v)[k+1] = ((*v)[k+1] - alpha) / r
		for j := k + 2; j < n; j++ {
			(*v)[j] /= r
		}

		S2OPAI(n, v, P)
		matMul(n, P, A, TEMP, TPOSE)
		matMul(n, TEMP, P, A, TPOSE)
	}
}

func householder(
	len int,
	A *[][]float64,
	A1 *[][]float64,
	Q *[][]float64,
	H *[][]float64,
	TEMP *[][]float64,
	TPOSE *[][]float64,
	v *[]float64,
) {
	copy2DSlice(A1, A)

	s := (*A1)[len-1][len-1]
	t := (*A1)[len-2][len-2]
	x := (*A1)[len-2][len-1]
	d := (t - s) / 2
	shift := s - x*x/(d+sgn(d)*math.Sqrt(d*d+x*x))

	for i := 0; i < len; i++ {
		(*A1)[i][i] -= shift
	}

	for k := 0; k < len-1; k++ {
		for j := 0; j < k; j++ {
			(*v)[j] = 0
		}
		for j := k; j < len; j++ {
			(*v)[j] = (*A1)[j][k]
		}

		(*v)[k] += norm2(len, v) * sgn((*v)[k+1])
		vDivS(len, v, norm2(len, v))

		if k == 0 {
			S2OPAI(len, v, Q)
			matMul(len, Q, A1, TEMP, TPOSE)
			A1, TEMP = TEMP, A1
		} else {
			S2OPAI(len, v, H)
			matMul(len, H, A1, TEMP, TPOSE)
			A1, TEMP = TEMP, A1
			matMul(len, Q, H, TEMP, TPOSE)
			Q, TEMP = TEMP, Q
		}
	}

	matMul(len, A1, Q, A, TPOSE)

	for i := 0; i < len; i++ {
		(*A)[i][i] += shift
	}
}

func eigQR(
	iter int,
	A *[][]float64,
	eigVals *[]float64,
	T1 *[][]float64,
	T2 *[][]float64,
	T3 *[][]float64,
	T4 *[][]float64,
	T5 *[][]float64,
	v *[]float64,
) {
	for ind := len(*A); ind > 1; ind-- {
		for i := 0; i < iter; i++ {
			householder(ind, A, T1, T2, T3, T4, T5, v)
		}
		(*eigVals)[ind-1] = (*A)[ind-1][ind-1]
	}

	(*eigVals)[0] = (*A)[0][0]
}

func eigenLoop(size int, iter int, loops int) {
	A := create2DSlice(size)

	T1 := create2DSlice(size)
	T2 := create2DSlice(size)
	T3 := create2DSlice(size)
	T4 := create2DSlice(size)
	T5 := create2DSlice(size)
	v := make([]float64, size)

	eigVals := make([]float64, size)

	for loop := 0; loop < loops; loop++ {
		randSymmetric(&A)
		tridiagonalize(size, &A, &T1, &T2, &T3, &v)
		eigQR(iter, &A, &eigVals, &T1, &T2, &T3, &T4, &T5, &v)
	}
}

func main() {
	if len(os.Args) != 5 {
		fmt.Println("Usage: ./eigenvalue_go <Matrix Size> <Convergence Iterations> <Loops> <Threads>")
		return
	}

	args := make([]int, 4)
	for i := 0; i < 4; i++ {
		arg, err := strconv.ParseInt(os.Args[i+1], 10, 64)
		if err != nil {
			panic(err)
		}
		args[i] = int(arg)
	}

	if args[3] == 1 {
		eigenLoop(args[0], args[1], args[2])
		os.Exit(0)
	}

	var wg sync.WaitGroup
	wg.Add(int(args[3]))

	for i := 0; i < args[3]; i++ {
		go func() {
			eigenLoop(args[0], args[1], args[2])
			wg.Done()
		}()
	}

	wg.Wait()
}
