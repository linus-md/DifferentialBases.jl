# DifferentialBases.jl

[![Build Status](https://github.com/linus-md/DifferentialBases.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/linus-md/DifferentialBases.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/github/linus-md/DifferentialBases.jl/graph/badge.svg?token=6XCVN0734M)](https://codecov.io/github/linus-md/DifferentialBases.jl)

DifferentialBases.jl is a Julia package for computing Groebner bases in differential settings.

### How to use DifferentialBases.jl

The following example is derived from a simple pendulum.

```julia
using DifferentialBases: differential_basis
using AlgebraicSolving: polynomial_ring, GF, Ideal
R, variables = polynomial_ring(GF(101),["dl","x","y","u","v",""], internal_ordering=:degrevlex)
(dl,x,y,u,v,l) = variables
derivatives = Dict(x => u, y => v, u => x*l, v => y*l - 1, l => dl)
ideal = Ideal([x^2 + y^2 - 1])
```

Calling `differential_basis(ideal, derivatives, R, 0)` results in the following Gröbner basis:

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
