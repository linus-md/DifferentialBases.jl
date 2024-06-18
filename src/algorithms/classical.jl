using AbstractAlgebra: vars
using AlgebraicSolving: polynomial_ring, Ideal, GF

function diff_op(f, derivatives)
    # TODO
end

function intersect(I, S)
    # This needs a generators in form of a Gröbner basis with an elimination order
    sub_ideal = []
    for generator in I.gens
        symbols = [Symbol(var) for var in vars(generator)]
        if issubset(symbols, S.data.S)
            push!(sub_ideal, generator)
        end
    end
    return sub_ideal
end

function differential_basis(ideal, derivatives)
    # TODO
    # Infer the subring
end

R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:degrevlex)
I = Ideal([x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y, x*y])

S, (x,y) = polynomial_ring(GF(101), ["x","y"], internal_ordering=:lex)
intersect(I, S)
