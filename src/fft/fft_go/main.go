package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"sync"
)

func fft(x *[]float64) {
	var real, imag, tReal, tImag float64

	n := len(*x)

	if n == 2 {
		return
	}

	even := make([]float64, n/2)
	odd := make([]float64, n/2)

	for i := 0; i < n/2; i += 2 {
		even[i] = (*x)[2*i]
		even[i+1] = (*x)[2*i+1]
		odd[i] = (*x)[2*i+2]
		odd[i+1] = (*x)[2*i+3]
	}

	fft(&even)
	fft(&odd)

	for k := 0; k < n/2; k += 2 {
		imag = -2 * math.Pi * float64(k) / float64(n)
		real = math.Cos(imag)
		imag = math.Sin(imag)

		tReal = real*odd[k] - imag*odd[k+1]
		tImag = real*odd[k+1] + imag*odd[k]

		(*x)[k] = even[k] + tReal
		(*x)[k+1] = even[k+1] + tImag

		(*x)[k+n/2] = even[k] - tReal
		(*x)[k+n/2+1] = even[k+1] - tImag
	}
}

func fftLoop(size int64, loops int64) {
	x := make([]float64, 2*size)

	for loop := int64(0); loop < loops; loop++ {
		for i := range x {
			x[i] = rand.Float64()
		}
		fft(&x)
	}
}

func main() {

	if len(os.Args) != 4 {
		fmt.Println("Usage: ./fft_go <Input Size> <Loops> <Threads>")
		return
	}

	args := make([]int64, 3)
	for i := 0; i < 3; i++ {
		arg, err := strconv.ParseInt(os.Args[i+1], 10, 64)
		if err != nil {
			panic(err)
		}
		args[i] = arg
	}

	args[0] = 1 << args[0]

	if args[2] == 1 {
		fftLoop(args[0], args[1])
		os.Exit(0)
	}

	var wg sync.WaitGroup
	wg.Add(int(args[2]))

	for i := int64(0); i < args[2]; i++ {
		go func() {
			fftLoop(args[0], args[1])
			wg.Done()
		}()
	}

	wg.Wait()

}
