# MFE1 Input File Format

The MFE1 input file is a [JSON](https://json.org) document that consists of four
primary sublists:
```
{
    "mfe-model": model-sublist,
    "initial-condition": ic-sublist,
    "simulation": simulation-sublist,
    "solver": solver-sublist
}
```
These sublists are described in the following.

## MFE Model Specification
The MFE model sublist contains the input parameter that define the PDE system
and MFE algorithm.

*model-sublist* is:
```
{
    "pde-library": string,
    "pde-params": pde-params-sublist,
    "boundary-conditions": bc-sublist,
    "kreg": integer,
    "eltvsc": float OR float-N-vector,
    "segspr": float OR float-N-vector,
    "fdinc": float,
    "eqw": float-N-vector // optional
}
```
**pde-library**

The name of the shared object library that implements the parts of the
discrete MFE equations that are specific to a particular PDE or PDE system.
The following are included with the MFE1 installation:
- "libconv-diff.so": convection-diffusion equation (N=1)
- "libburgers.so": viscous Burgers' equation (N=1)
- "libnavier-stokes.so": compressible Navier-Stokes equations (N=3)
- "libdrift-diff.so": single-carrier semiconductor drift-diffusion equation
  (N=2)

Several of the input parameters, here and in other sections, either take a
single scalar value in the case of a scalar PDE (N=1) or a vector value
of length N equal to the number of PDEs in the case of a system of PDE (N>1).

**kreg**  
Specifies the type of dynamic regularization to use:
* 1 for rate of deformation regularization; and
* 2 for total gradient regularization.

See section 2.6 of [1] for details.

**eltvsc**  
The viscosity coefficient ($A^2$) in the dynamic regularization; see section
2.6 of [1].

**segspr**  
The compressive spring force coefficient ($C^2$) in the static regularization;
see section 2.7 of [1].

**fdinc**  
The increment used when computing the finite difference approximation to the
Jacobian of the discrete MFE system.

**eqw**  
The equation weights in the discrete MFE system. This is only relevant to PDE
systems (N>1), and defaults to equal weighting (1) of the equations in the
system.

### PDE Parameters
The pde-params-sublist contains the input parameters that are specific to a
particular PDE system. It is consumed by the source code for the PDE library;
refer to that code for details on what the sublist is expected to define.

For the included libraries, *pde-params-sublist* is:
* for "libconv-diff.so", "libburgers.so", "libnavier-stokes.so"
  ```
  { "visc": float }
  ```
  where **visc** is the diffusion coefficient.

* for "libdrift-diff.so"
  ```
  {
    "lambda": float, "eps": float,
    "u-scale-factor": float, "v-scale-factor": float
  }

### Boundary Condition Specification
The options for boundary conditions are extremely limited, either constant
Dirichlet, with boundary value drawn from the initial conditions, or
homogeneous Neumann.

*bc-sublist* is:
```
{
    "dir-left": boolean OR boolean-N-vector,
    "dir-right": boolean OR boolean-N-vector,
    "fixed-node-left": boolean, // optional, default true
    "fixed-node-right": boolean // optional, default true
}
```

**dir-left**, **dir-right**  
Specify whether a Dirichlet condition is imposed (true) on the corresponding
unknown at the respective end point, or a homogeneous Neumann condition is
imposed (false).

**fixed-node-left**, **fixed-node-right**  
By default the end point nodes are fixed (true), but they can be allowed to
move with velocities determined by the PDE system if assigned the value false.

## Initial Condition Specification
The initial continuous piecewise-linear finite element discretization of the
PDE unknowns is linearly interpolated from a "coarse grid" table of values
specified by this sublist.

*ic-sublist* is:
```
{
    "useg": float-matrix,
    "niseg": integer-vector,
    "ratio": float-vector // optional
}
```
The table of values is specified by **useg**. It takes the form
```
[ float-vector, float-vector, ... ]
```
where each *float-vector* is of the form [ $x$, $u_1$, ..., $u_N$]
and specifies a coarse grid "point". The intervals between successive points
are subdivided into subintervals, the number of which are specified by
**niseg**. By default the subintervals are equal sized, however if **ratio**
is specified, its values give the size ratio of successive subintervals within
each interval. In either case, the finite element discretization of the PDE
unknowns is linearly interpolated as functions of $x$ between their coarse grid
values.

## Simulation Specification
This simulation sublist contains the parameters that control the simulation.

*simulation-sublist* is:
```
{
    "t-init": float,
    "h-init": float,
    "output-times": float-vector,
    "max-step": integer, // optional
    "output-freq": integer, // optional
}
```
<!---"hard-event-lookahead": integer, // optional --->

**t-init**  
The initial problem time.

**h-init**  
The initial time step size.

**output-times**  
A list of problem times when the solution is to be output. The listed values
need not be ordered. Moreover any that do not exceed the initial problem time
are ignored. Note that the time stepping will not hit these times exactly.
Instead, a high-order interpolated solution, using the solution at recent time
steps, is output at the first time step that reaches or surpasses the requested
output time. The simulation is finished when the last output time (i.e., the
greatest) is reached.

**max-step**  
The maximum number of time steps allowed. The simulation is halted if the
number of time steps exceeds this value before the last output time is
reached. The default is unlimited time steps.

**output-freq**  
Specifies how frequently, in number of steps, the solution is output. This
is in addition to the output specified by **output-times**. The default is 0,
meaning no output frequency.

## Solver Specification
This sublist contains the parameters associated with the DAE integrator,
which uses the variable step-size BDF2 method. Each implicit step of BDF2
requires the solution of a nonlinear system and the nonlinear Krylov (NLK)
iterative method is used to solve it; see [1, §9] and [2].

The iterative method is deemed converged when the norm of the solution
correction $(\delta x, \delta u_1, \ldots, \delta u_N)$ is sufficiently
small, namely

$$ \max\left\\{
\frac{\lVert\delta x\rVert_\infty}{\epsilon_x},
\frac{\lVert\delta u_1\rVert_\infty}{\epsilon_1}, \ldots,
\frac{\lVert\delta u_N\rVert_\infty}{\epsilon_N},
\frac{\lVert\delta\Delta x/\Delta x\rVert_\infty}{\epsilon_r}
\right\\} < \eta, $$

where $\Delta x$ are the finite element cell sizes.

*solver-sublist* is:
```
{
    "abs-u-tol": float-N-vector,
    "abs-x-tol": float,
    "rel-dx-tol": float,
    "h-min": float,              // optional
    "h-max": float,              // optional
    "nlk-max-iter": integer,     // optional
    "nlk-tol": float,            // optional
    "nlk-max-vec": integer,      // optional
    "nlk-vec-tol": float,        // optional
    "max-try": integer,          // optional
    "verbose-stepping": boolean, // optional
    "dxmin": float               // optional
}
```

**abs-u-tol**, **abs-x-tol**, **rel-dx-tol**  
These specify the error tolerances $\epsilon_i$, $\epsilon_x$, and $\epsilon_r$
respectively.

**nlk-max-iter**  
The maximum number of iterations allowed for a nonlinear solve. If the iteration
does not converge within this number of iterations, the DAE integrator will
attempt to recover by repeating the step with an updated preconditioner or by
reducing the time step size. Failures of this sort are entirely normal and not
fatal unless repeated excessively. The default is 5 iterations.

**nlk-tol**  
Specifies the error tolerance $\eta$. The solution correction norm is defined
(via the $\epsilon$ tolerances) such that 1 is the target value for the
predictor error which drives the choice of time step size. We only need to
solve the nonlinear system to somewhat greater accuracy, thus the default
value is 0.1.

**nlk-max-vec**  
The maximum number of vectors used by the nonlinear Krylov acceleration (NKA)
algorithm. The default is to use all the solution iterates, namely one less
than the value of **nlk-max-iter**. A value of zero will disable use of NKA
and the resulting iteration is reduced to a (preconditioned) fixed point
iteration.

**nlk-vec-tol**  
The vector drop tolerance in the NKA algorithm. The default is 0.01.

**h-min**  
The minimum allowed time step size. Time stepping and the simulation are
halted if the DAE integrator wants to use a smaller step size. The default
is the smallest representable 64-bit floating point number. It is a good idea
to set this variable or **max-step** to avoid run-away simulations. 

**h-max**  
The maximum allowed time step size. The DAE integrator will never use a
larger step size even if it would otherwise want to do so. The default
is no limit.

**max-try**  
The maximum number of attempts at a single time step. When a time step fails
for any reason, the DAE integrator will try to recover, usually by reducing
the step size. If it does not succeed within this number of attempts,
time stepping and the simulation are halted. The default value is 10.

**verbose-stepping**  
If set to true, copious diagnostic information about the time stepping
procedure is written to the file `bdfout`.

**dxmin**  
The DAE integrator checks that all solution iterates are *admissable*. For MFE1
this means that all cell sizes are no smaller than the lower bound specified
by this variable. Its default value is the smallest representable 64-bit
floating point number. 

----

## References

1. N. N. Carlson and K. Miller. Design and Application of a Gradient-Weighted
    Moving Finite Element Code I: in One Dimension. SIAM Journal on Scientific
    Computing, 19(3):728-765, 1998.

2. Nonlinear Krylov Acceleration, https://github.com/nncarlson/nka
