# Compiler Efficiency Suite
### An Evaluation of the Energy Characteristics of General-Purpose Programming Languages

This repository contains the benchmark suite and build system developed for my undergraduate thesis.
The benchmark suite consists of three benchmarks, each of which was translated into the following lanaguages:  
`C`, `C#`, `C++`, `Go`, `Java`, `Python`, `Ruby`, `Rust`, and `Swift`

## Benchmarks
1. #### Eigenvalue
   CLI Usage: `./eigenvalue <Matrix Size> <Convergence Iterations> <Loops> <Threads>`
    * `<Matrix Size>` - Side length (N) of the input matrices (N x N).
    * `<Convergence Iterations>` - Number of iterations toward convergence per eigenvalue.
    * `<Loops>` - Number of input matrices to run per thread.
    * `<Threads>` - Number of concurrent threads.
2. #### FFT
   CLI Usage: `./fft <Input Size> <Loops> <Threads>`
    * `<Input Size>` - Length of the complex input arrays (2^N).
    * `<Loops>` - Number of input arrays to run per thread.
    * `<Threads>` - Number of concurrent threads.
3. #### SHA512
   CLI Usage: `./sha512 <File Path> <Iterations>`
    * `<File Path>` - Path to target file for hashing.
    * `<Iterations>` - Number of SHA512 block hashing iterations to perform per input block.

## Requirements
The benchmark suite and build system were designed with the following assumptions:
- Linux installation with ring 0 access
- Access to the power/energy-pkg/ counter through perf
- An sh-compatible shell like Bash

## Building
The build system designed for this benchmark suite is configured through the config.json file.
Before building, compiler binaries must be in the paths specified in the config.json file.
This repository hosts an example config.json with the same configuration options that were used 
to evaluate compiler efficiencies in my thesis.
To build the suite, execute the following commands:
```
python not_quite_autogen.py
./build.sh
```
Additionally, you may want to create files to hash in *sha512*.
This can be done using the file generator tool in the `tools/filegen/` folder.
To create the files used in the thesis, use the following commands after compiling filegen:
```
./filegen test.txt 10000000000 False
./filegen test_small.txt 30000000 False
```
*Note: This will create a 10 GB file on your drive. You have been warned.*

## Running
With root privileges, run: 
```
./benchmark.sh
```
This will generate a results folder containing time and energy values for each benchmark run.

## Parsing Results
From the project root, run:
```
python tools/resultsparser/parse.py
```  
This will generate a results_parsed folder containing csv files for each benchmark configuration tested.  
*Note: The alternate benchmark configurations in the parse.py file do not sync with the build system configuration, 
and they need to be updated accordingly if they are modified.*