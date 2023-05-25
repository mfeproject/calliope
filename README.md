# Calliope

**This is a work in progress.**  My immediate goal is to create general 1D, 2D,
and 3D GWMFE reference codes, implemented in modern object-oriented Fortran
that will replace the original archaic Fortran 77 codes (1D and 2D
[here](https://github.com/mfeproject/legacy-gwmfe)). These will be derived
from my existing unreleased Fortran 90 codes (1D is now available
[here](https://github.com/mfeproject/mfe1)).

## Status (May 2023)
The modernization of the 1D code is essentially complete and the code has
feature-parity with the legacy GWMFE1DS code. Some rough edges still remain
to be smoothed out over time. Development focus is now on enhancing
capabilities, algorithm improvements, and coarray parallelization.

### Previous Updates
**April 2023.**
After a nearly 8-year hiatus, I've returned to this in earnest. I'm currently
working on restructuring the 1D code with an eye to parallelizing the code
using coarrays (rather than openMP, for example) for multicore shared memory
systems. (1D is ridiculously fast in serial as it is, but I'm anticipating
this as a pattern for the 2D and 3D codes.)

## Compiling
The code depends on the petaca library whose source is available
 [here](https://github.com/nncarlson/petaca). After compiling and installing the
 library, set your environment variable `PETACA_ROOT` to your installation
 path before proceeding.

```shell
$ git clone https://github.com/mfeproject/calliope.git
$ cd calliope
$ mkdir build
$ cd build
$ cmake .. -D CMAKE_BUILD_TYPE=Release
$ make
```

You may need to set your `FC` environment variable to the path of your Fortran
compiler before running `cmake` to ensure that CMake finds the correct compiler.
The code works with the Intel oneAPI classic `ifort` compiler, GFortran 12.2,
and NAG 7.1. (GFortran 13.1 is buggy; *do not use.*)

This produces a single executable `mfe1` and several shared libraries (`*.so`
files) located in `build/src`, which can then be installed:

```shell
$ make install
```
The default installation location is `/opt/calliope`. This can be changed by
specifying a different location with a `-D CMAKE_INSTALL_PREFIX=/my/location`
argument to the above `cmake` command.

## Testing
From the `build` directory give the command `ctest` to run the tests. The tests
run the code and compare against reference results.

## Examples
The `examples` directory contains input files for some example problems. To
run them you'll need to add `/opt/calliope/bin` to your PATH environment
variable and invoke the `mfe1` program with an input file as an argument;
for example
```shell
mfe1 sod-shock-tube.json
```
When `mfe1` is executed it will load one of the shared libraries that provides
procedures specific to a particular system of PDEs. The current PDE systems are:

* Scalar convection-diffusion equation
* Burgers' equation
* Single-carrier semiconductor drift-diffusion equations
* Compressible Navier-Stokes equations

The program writes diagnostic output to `mfelog` and the solution at a sequence
of time steps to `mfegrf`. The format of this file is fairly self-explanatory.
Many tools are capable of visualizing this 1D data. But note that the custom
graphics program `gp1` from the legacy GWMFE1DS code,
https://github.com/mfeproject/legacy-gwmfe/tree/master/gwmfe1ds,
is able to read this file.

## Documentation
See [doc/input-format.md](doc/input-format.md) for documentation of the input
file format.
