using AbstractAlgebra: vars, derivative
using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis, normal_form

"""
    partial(q, derivatives)

    This function evaluates the linear differential operator from 
    Definition 7 and Notation 3 in the thesis for a polynomial and derivatives.

    # Arguments
    - `q`: a polynomial
    - `derivatives`: a dictionary of derivatives

    # Returns
    - the linear differential operator applied to the polynomial
"""
function partial(q, derivatives)    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n "There can't be more derivatives than variables."
    
    result = 0
    for (var, value) in derivatives
        result += value * derivative(q, var)
    end
    return result
end

"""
    cap(G, S_vars)

    This function computes the intersection of a Groebner basis with a set of variables.

    # Arguments
    - `G`: a Groebner basis
    - `S_vars`: a set of variables

    # Returns
    - the intersection of the Groebner basis with the set of variables
"""
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


"""
    differential_basis(ideal, derivatives, R, nf=false, info_level=0)

    This function computes the differential Gröbner basis of an ideal or 
    an approximation of a full differential ideal if not all derivatives
    are known.

    # Arguments
    - `ideal`: an ideal
    - `derivatives`: a dictionary of derivatives
    - `R`: a polynomial ring
    - `nf`: a boolean indicating whether to compute the normal form
    - `info_level`: an integer indicating the level of information to print

    # Returns
    - the differential Gröbner basis
"""
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
