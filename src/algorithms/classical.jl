using AbstractAlgebra: vars, derivative
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis, normal_form

function partial(q, derivatives)
    # The i-th derivative must correspond to the i-th variable, 
    # all trailing variables are not derived.
    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n "There can't be more derivatives than variables."
    
    result = 0
    # TODO use correct derivative
    for (var, value) in derivatives
        result += value * derivative(q, var)
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
    S_vars = [var for (var, _) in derivatives]
    i = 1
    println("i = ", i)

    # Start computing the differential basis
    G1 = groebner_basis(I)
    G1_S = intersect_ideal(G1, S_vars)
    G1_S = [partial(g, derivatives) for g in G1_S]
    append!(G1, G1_S)
    G2 = groebner_basis(Ideal(G1))

    # Repeat until closed under partial
    while Ideal(G1) != Ideal(G2)
        i += 1
        println("i = ", i)
        G1 = G2
        G1_S = intersect_ideal(G1, S_vars)
        G1_S = [partial(g, derivatives) for g in G1_S]
        append!(G1, G1_S)
        G2 = groebner_basis(Ideal(G1))
    end
    return G2
end

# invlex order
R, (dl,l,v,u,y,x) = polynomial_ring(GF(101),["dl","l","v","u","y","x"], internal_ordering=:lex)
derivatives = Dict(
    x => u,
    y => v,
    u => x*l,
    v => y*l - 1,
    l => dl
)

I = Ideal([x^2 + y^2 - 1])
#differential_basis(I, derivatives)
