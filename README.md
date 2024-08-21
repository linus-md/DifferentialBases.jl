# DifferentialBases.jl

[![Build Status](https://github.com/linus-md/DifferentialBases.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/linus-md/DifferentialBases.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/github/linus-md/DifferentialBases.jl/graph/badge.svg?token=6XCVN0734M)](https://codecov.io/github/linus-md/DifferentialBases.jl)

DifferentialBases.jl is a Julia package for computing Groebner bases in differential settings.

### How to use DifferentialBases.jl

The following example is derived from a simple pendulum. Calling ``simple_pendulum()`` loads the relevant equations.

```julia
using DifferentialBases: differential_basis, simple_pendulum
using AlgebraicSolving: polynomial_ring, GF, Ideal
ideal, derivatives, R =  simple_pendulum()
```

Then we have the following system:

```julia-repl
println(ideal)
Nemo.FqMPolyRingElem[x^2 + y^2 + 100]
println(derivatives)
Dict{Nemo.FqMPolyRingElem, Nemo.FqMPolyRingElem}(y => v, x => u, u => x*l, v => y*l + 100, l => dl), Multivariate polynomial ring in 6 variables over GF(101)
```

Calling `differential_basis(ideal, derivatives, R, false, 0)` results in the following Gröbner basis:

```julia
8-element Vector{Nemo.FqMPolyRingElem}:
 y^2 + x^2 + 100
 v*y + u*x
 l + v^2 + u^2 + 100*y
 l*x^2 + 100*l + 100*u^2 + 100*y*x^2 + y
 100*v*x^2 + v + u*y*x
 l*y*x + 100*v*u + x^3 + 100*x
 l*y + 100*v*u*x + u^2*y + x^2 + 100
 dl + 98*v
```

We can uncover interesting additional constraints from those equations i.e $ux + vy = 0$ describes that the motion of the masspoint of the pendulum is tangent to the rod of the pendulum. Furthermore we uncover that $u^2 + v^2 - y + 1 = 0$ which describes that the sum of potential and kinetic energy is constant i.e. energy is conserved. We learn that $dl= 3v$.

### Installation

To install DifferentialBases.jl enter the **Pkg** REPL by pressing `]` and execute

```julia-repl
add https://github.com/linus-md/DifferentialBases.jl
```
