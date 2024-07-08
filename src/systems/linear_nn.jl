using AlgebraicSolving: polynomial_ring, GF, Ideal
using AbstractAlgebra: derivative, matrix

function linear_nn_2(m, n, r)
    R, A, B, x, y = polynomial_ring(
        GF(101),
        :A=>(1:m, 1:n),
        :B=>(1:n, 1:r),
        :x => (1:r),
        :y => (1:m))

    nvars = m*n + n*r + r + m
    A, B = matrix(A), matrix(B)

    f = A*B*x - y
    f = f[1]^2 + f[2]^2
    
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
    return ideal, derivatives, R
end

function linear_nn_3(m, n, r, s)
    R, A, B, C, x, y = polynomial_ring(
        GF(101),
        :A=>(1:m, 1:n), 
        :B=>(1:n, 1:r),
        :C=>(1:r, 1:s),
        :x => (1:s), 
        :y => (1:m))

    nvars = m*n + n*r + r*s + s+m
    A, B, C = matrix(A), matrix(B), matrix(C)

    f = A*B*C*x - y
    f = f[1]^2 + f[2]^2
    
    derivatives = Dict(
        R[i] => derivative(f, i) for i in 1:nvars - s - m + 1)
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
    ideal = Ideal(D_flat)
    return ideal, derivatives, R
end
