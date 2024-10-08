using AlgebraicSolving: polynomial_ring, QQ, Ideal, groebner_basis, normal_form
using AbstractAlgebra: derivative, matrix, vars

function activation_function(x)
    return x^2
end


"""
    poly_nn_2(m, n, r)

Implementation of example 2.5.9 in the thesis.

# Arguments
- `m`: number of rows in A
- `n`: number of columns in A and rows in B
- `r`: number of columns in B

Returns the corresponding system of the two layer polynomial neural network
"""
function poly_nn_2(m, n, r, activation=activation_function)
    R, A, B, x, y = polynomial_ring(
        QQ,
        :A => (1:m, 1:n),
        :B => (1:n, 1:r),
        :x => (1:r),
        :y => (1:m))

    nvars = m*n + n*r + r + m
    A, B = matrix(A), matrix(B)

    f = A*activation.(B*x) - y
    f = sum(f.*f)

    derivatives = Dict(
        R[i] => derivative(f, i) for i in 1:nvars - r - m + 1)
    for i in nvars - r - m + 1:nvars
        derivatives[R[i]] = R(0)
    end
    
    C = transpose(A) * A - B *  transpose(B)
    C_flat = eltype(A)[]
    for elem in C
        push!(C_flat, elem)
    end
    ideal = Ideal(C_flat)
    return ideal, derivatives, R, []
end
