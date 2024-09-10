using AlgebraicSolving
using AbstractAlgebra

"""
    linear_nn_2(m, n, r)

Implementation of example 2.5.7 in the thesis.

# Arguments
- `m`: number of rows in A
- `n`: number of columns in A and rows in B
- `r`: number of columns in B

Returns the corresponding system of the two layer linear neural network
"""
function linear_nn_2(m, n, r)    
    R, A, B, x, y = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.QQ,
        :A=>(1:m, 1:n),
        :B=>(1:n, 1:r),
        :x => (1:r),
        :y => (1:m))

    nvars = m*n + n*r + r + m
    A, B = AbstractAlgebra.matrix(A), AbstractAlgebra.matrix(B)

    f = A*B*x - y
    f = sum(f.*f)
    
    derivatives = Dict(
        R[i] => AlgebraicSolving.derivative(f, i) for i in 1:nvars - r - m + 1)
    for i in nvars - r - m + 1:nvars
        derivatives[R[i]] = R(0)
    end
    
    C = transpose(A) * A - B *  transpose(B)
    C_flat = eltype(A)[]
    for elem in C
        push!(C_flat, elem)
    end
    ideal = AlgebraicSolving.Ideal(C_flat)
    return ideal, derivatives, R, []
end

"""
    linear_nn_3(m, n, r, s)

Implementation of example 2.5.8 in the thesis.

# Arguments
- `m`: number of rows in A
- `n`: number of columns in A and rows in B
- `r`: number of columns in B and columns in C
- `s`: number of rows in C

Returns the corresponding system of the three layer linear neural network
"""
function linear_nn_3(m, n, r, s)
    R, A, B, C, x, y = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.QQ,
        :A=>(1:m, 1:n), 
        :B=>(1:n, 1:r),
        :C=>(1:r, 1:s),
        :x => (1:s), 
        :y => (1:m))

    nvars = m*n + n*r + r*s + s+m
    A = AbstractAlgebra.matrix(A)
    B = AbstractAlgebra.matrix(B)
    C = AbstractAlgebra.matrix(C)

    f = A*B*C*x - y
    f = sum(f.*f)
    
    derivatives = Dict(
        R[i] => AlgebraicSolving.derivative(f, i) for i in 1:nvars - s - m + 1)
    for i in nvars - s - m + 1:nvars
        derivatives[R[i]] = R(0)
    end
    
    D1 = transpose(A) * A - B *  transpose(B)
    D2 = transpose(B) * B - C *  transpose(C)
    D_flat = eltype(A)[]
    for elem in D1
        push!(D_flat, elem)
    end
    for elem in D2
        push!(D_flat, elem)
    end
    ideal = AlgebraicSolving.Ideal(D_flat)
    return ideal, derivatives, R, []
end
