# Branch `pdt`

This was an interesting and nontrivial experiment to see if the `NodeVar` and
`NodeMtx` derived types could be parametrized by the number of PDEs (NEQN)
rather than hardwired with that value.

This actually works! But only with the Intel oneAPI 2021.9 compiler (ifort).
The NAG compiler fails fatally at two points (at least); see
[nag-20230426.f90](https://github.com/nncarlson/fortran-compiler-tests/blob/master/nag-bugs/nag-20230426.f90) and
[nag-20230426b.f90](https://github.com/nncarlson/fortran-compiler-tests/blob/master/nag-bugs/nag-20230426b.f90).
GFortran was expected to fail, and did, due to its completely botched PDT
implementation in 12.2 and earlier. There has been some very recent work on
PDTs so the situation may improve soon.

The `mfe_types.f90` defined a collection of binary operators for the types,
however only multiplication by a constant and assignment of a scalar real were
being used, and only those were modified; the rest were commented out. To start
these were left as generics that overloaded * and =.

`time` for the Sod test problem using the original code reports 0.10 sec

As anticipated, this PDT version is slower, but at 0.98 sec the amount of
slowdown is surprising.

I then made the binary operations type bound rather than "loose". Fortunately
this works as well and is essentially the same speed at 0.99 sec.

It will be interesting to see if NAG (once fixed) shows the same extreme
slow down that Intel shows.

Regardless, this data organization -- array-of-structs -- should be replaced
by a struct-of-arrays or even plain rank-1 or rank-2 array organization,
which are generally considered to be more performant.
