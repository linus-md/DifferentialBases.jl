# DifferentialBases.jl

[![Build Status](https://github.com/linus-md/DifferentialBases.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/linus-md/DifferentialBases.jl/actions/workflows/CI.yml?query=branch%3Amain)

DifferentialBases.jl is a Julia package for computing Groebner bases in differential settings. 

### How to use DifferentialBases.jl

The following example is derived from a simple pendulum. 

```julia
using AbstractAlgebra
R, (x1, x2, x3, x4, x5) = QQ["x1", "x2", "x3", "x4", "x5"]
ideal = [x1^2 + x2^2 - 1]
partial = [x3, x4, x5*x1, x5*x2 - 1]
```

Calling `differential_bais(ideal, partial)` results in the following Gröbner basis:
