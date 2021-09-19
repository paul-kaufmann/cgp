## About

This is an adaptation of Julian Miller's Cartesian Genetic Programming (CGP) implementation with new Evolutionary Algorithms and Benchmarks.

## Requirements

* g++-11

## How to Use

There are two independent code trees for Boolean and regression benchmarks. In each implementation there is a Makefile in the source-files subdirectory. Calling it will compile the cgp and cgp-boolean executables. The executables will show a complete list of parameters. The parameter files are stored in the parameters directory. For Boolean benchmarks one need a set of two files defining the arithmetic functions (eg and-or-nand-nor.par) and the input-output mappings (e.g. add3.plu). For regression benchmarks one need instead of a single input-output mappings definition a definition for training and a definition for test samples.
