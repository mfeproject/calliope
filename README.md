# Calliope

**This is a work in progress.**  My immediate goal is to create general 1D, 2D,
and 3D GWMFE reference codes, implemented in modern object-oriented Fortran
that will replace the original archaic Fortran 77 codes (1D and 2D
[here](https://github.com/mfeproject/legacy-gwmfe)). These will be derived
from my existing unreleased Fortran 90 codes (1D is now available
[here](https://github.com/mfeproject/mfe1)).

### Status (April 2023)
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
$ cmake .. -DCMAKE_BUILD_TYPE=Release
$ make
```

You may need to set your `FC` environment variable to the path of your Fortran
compiler before running `cmake` to ensure that CMake finds the correct compiler.
The code works with the Intel oneAPI classic `ifort` compiler, GFortran 12.2,
and NAG 7.1.

For the purposes of development there is a single executable `gas` located in
`build/src` that solves the 1D Navier-Stokes equations.

## Testing
From the `build` directory give the command `ctest` to run the tests. The only
test at present is the Sod shock tube problem, which runs the code and compares
against reference results.

## Examples
The `test/sod` directory contains an input file `mfein` for the Sod shock tube
problem; see [1].  Assuming the executable `gas` is in your path,
```shell
gas mfein
```
will run the problem. Diagnostic output is written to `mfelog` and the solution
at a sequence of time steps is written to `mfegrf`. The format of this file is
fairly self-explanatory. The graphics program `gp1` from the legacy GWMFE1DS
code, https://github.com/mfeproject/legacy-gwmfe/tree/master/gwmfe1ds, is able
to read this file.

---

[1] N. N. Carlson and K. Miller, "Design and Application of a Gradient-Weighted Moving Finite Element Code I: in One Dimension", SIAM J. Sci. Comput., 19(3), 728-765 (1998).

