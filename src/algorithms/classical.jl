using AbstractAlgebra: vars, derivative
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis, normal_form

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

function intersect_ideal(G, S_vars)
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

function differential_basis(I, derivatives)
    # Infer and create the subring
    k = length(derivatives)
    S_vars = R.data.S[1:k]
    i = 1
    println(i)

    # Start computing the differential basis
    G1 = groebner_basis(I)
    G1_S = intersect_ideal(G1, S_vars)
    G1_S = [normal_form(partial(g, derivatives), Ideal(G1)) for g in G1_S]
    G2 = groebner_basis(Ideal(append!(G1, G1_S)))

    while G1 != G2
        i += 1
        println(i)
        G1 = G2
        G1_S = intersect_ideal(G1, S_vars)
        G1_S = [normal_form(partial(g, derivatives), Ideal(G1)) for g in G1_S]
        G2 = groebner_basis(Ideal(append!(G1, G1_S)))
    end
    return G2
end

R, (dl,x,y,u,v,l) = polynomial_ring(GF(101),["dl","x","y","u","v","l"], internal_ordering=:lex)
derivatives = [0, u, v, x*l, y*l -1, dl]
I = Ideal([x^2 + y^2 - 1])
differential_basis(I, derivatives)
