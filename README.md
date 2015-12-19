## Work in Progress
This is a work in progress.  I have two near-term goals.  The first is to
create general 1D, 2D, and 3D GWMFE reference codes, implemented in modern
object-oriented Fortran that will replace the original archaic Fortran 77
codes.  These will be derived from my existing unreleased Fortran 90 codes.
The second is to have a suitable code base where I can experiment with
different vectorization and threading strategies for many-core architectures
like the Intel Xeon Phi.

#### Status (Dec 2015)
I'm working on "objectizing" the BDF2 integrator, and have a simple 1D
Galerkin heat equation code as a test driver. The repository includes Petaca
as a submodule, which is where I've been developing some basic infrastructure
components over the last few years. The next step is to pull in all the
specific code to generate the semi-discrete GWMFE nonlinear ODE system (in 1D).

## Cloning the repository
After cloning the repository you'll need to fetch the petaca submodule:

    $ git clone https://github.com/mfeproject/calliope.git
    $ cd calliope
    $ git submodule update --init
