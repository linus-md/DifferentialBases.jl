using AbstractAlgebra: vars, derivative
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis, normal_form, QQ, eliminate

function partial(q, derivatives)    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n "There can't be more derivatives than variables."
    
    result = 0
    # TODO use correct derivative
    for (var, value) in derivatives
        result += value * derivative(q, var)
    end
    return result
end

function cap(G, S_vars)
    # This needs a Gröbner basis with an elimination order
    sub_ideal = []
    for generator in G
        symbols = [Symbol(var) for var in vars(generator)]
        if issubset(symbols, S_vars)
            push!(sub_ideal, generator)
        end
    end
    return sub_ideal
end

function differential_basis(ideal, derivatives, R, verbose=false)
    # Infer and create the subring
    n = R.data.nvars
    S_vars = [Symbol(var.first) for var in derivatives]
    k = length(S_vars)
    eliminate = n - k
    if verbose
        i = 1
        println("i = ", i)
    end
    
    # Start computing the differential basis
    G1 = groebner_basis(ideal)
    pG1 = [partial(g, derivatives) for g in cap(G1, S_vars)]
    append!(pG1, G1)
    G2 = groebner_basis(Ideal(pG1), eliminate=eliminate)

    # Repeat until closed under partial
    while G1 != G2
        if verbose
            i += 1
            println("i = ", i)
        end
        G1 = G2
        pG1 = [partial(g, derivatives) for g in cap(G1, S_vars)]
        append!(pG1, G1)
        G2 = groebner_basis(Ideal(pG1), eliminate=eliminate)
    end
    return G1
end
