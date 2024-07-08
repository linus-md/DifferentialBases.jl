using AbstractAlgebra: vars, derivative
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis, normal_form

function partial(q, derivatives)    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n "There can't be more derivatives than variables."
    
    result = 0
    for (var, value) in derivatives
        result += value * derivative(q, var)
    end
    return result
end

function cap(G, S_vars)
    sub_ideal = []
    for generator in G
        symbols = [Symbol(var) for var in vars(generator)]
        if issubset(symbols, S_vars)
            push!(sub_ideal, generator)
        end
    end
    return sub_ideal
end

function differential_basis(ideal, derivatives, R, nf=false, info_level=0)
    # Infer and create the subring
    n = R.data.nvars
    S_vars = [Symbol(var.first) for var in derivatives]
    k = length(S_vars)
    eliminate = n - k
    
    # Start computing the differential basis
    G1 = groebner_basis(ideal)
    pG1 = [partial(g, derivatives) for g in cap(G1, S_vars)]
    if nf == true
        pG1 = [normal_form(pg, Ideal(G1)) for pg in pG1]
    end
    append!(pG1, G1)
    G2 = groebner_basis(Ideal(pG1), eliminate=eliminate,
                        intersect=false, info_level=info_level)
    if info_level > 0
        i = 1
        println("iteration ", i)
        println("#G = ", length(G1))
    end

    # Repeat until closed under partial
    while G1 != G2
        if info_level > 0
            i += 1
            println("iteration ", i)
            println("#G = ", length(G2))
        end
        G1 = G2
        pG1 = [partial(g, derivatives) for g in cap(G1, S_vars)]
        if nf == true
            pG1 = [normal_form(pg, Ideal(G1)) for pg in pG1]
        end
            append!(pG1, G1)
        G2 = groebner_basis(Ideal(pG1), eliminate=eliminate,
                            intersect=false, info_level=info_level)
    end
    return G1
end
