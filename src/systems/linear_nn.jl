using AlgebraicSolving: polynomial_ring, GF, Ideal
using AbstractAlgebra: derivative
#using LinearAlgebra

function linear_nn_2(m, n, r)
    R, A, B, x, y = polynomial_ring(
        GF(101),
        :A=>(1:m, 1:n),
        :B=>(1:n, 1:r),
        :x => (1:n),
        :y => (1:n))

    nvars = m*n + n*r + 2*n

    f = A*B*x - y
    f = f[1]^2 + f[2]^2
    
    derivatives = Dict(
        i => derivative(f, i) for i in 1:nvars - 2*n + 1)
    for i in nvars - 2*n + 1:nvars
        derivatives[i] = R(0)
    end
    
    #constraints = collect(Iterators.flatten(A*A' - B'*B))
    ideal = Ideal([R(0)])
    return ideal, derivatives
end

function linear_nn_3(m, n, r, s)
    R, A, B, C, x, y = polynomial_ring(
        GF(101),
        :A=>(1:m, 1:n), 
        :B=>(1:n, 1:r),
        :C=>(1:r, 1:s),
        :x => (1:n), 
        :y => (1:n))

    nvars = m*n + n*r + r*s + 2*n

    f = A*B*C*x - y
    f = f[1]^2 + f[2]^2
    
    derivatives = Dict(
        i => derivative(f, i) for i in 1:nvars - 2*n + 1)
    for i in nvars - 2*n + 1:nvars
        derivatives[i] = R(0)
    end
    
    #constraints = collect(Iterators.flatten(A*A' - B'*B))
    #push!(constraints, collect(Iterators.flatten(B*B' - C'*C)))
    ideal = Ideal(R(0))
    return ideal, derivatives
end

f = linear_nn_2(2, 2, 2)