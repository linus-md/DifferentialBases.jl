using AbstractAlgebra: vars
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis

function diff_op(f, derivatives)
    # TODO
end

function intersect(G, S)
    # This needs a generators in form of a Gröbner basis with an elimination order
    sub_ideal = []
    for generator in G
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
