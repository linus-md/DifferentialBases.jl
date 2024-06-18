using AbstractAlgebra: vars, derivative
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis

function partial(q, derivatives)
    # The i-th derivative must correspond to the i-th variable, 
    # all trailing variables are not derived.
    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n "There can't be more derivatives than variables."
    
    result = 0
    for (i, der) in enumerate(derivatives)
        result += der * derivative(q, i)
    end
    return result
end

function intersect_ideal(G, S)
    # This needs a Gröbner basis with an elimination order
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
