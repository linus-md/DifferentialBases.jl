using AbstractAlgebra
using AlgebraicSolving

"""
    diff_op(q, derivatives)

    This function evaluates the linear differential operator from 
    Definition 7 and Notation 3 in the thesis for a polynomial and derivatives.
    The inputs should not explicitly depend on time.

    # Arguments
    - `q`: a polynomial
    - `derivatives`: a dictionary of derivatives

    # Returns
    - the linear differential operator applied to the polynomial
"""
function diff_op(q, derivatives)    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n "There can't be more derivatives than variables."
    
    result = 0
    for (var, value) in derivatives
        result += value * AbstractAlgebra.derivative(q, var)
    end
    return result
end

"""
    intersect(G, S_vars)

    This function computes the intersection of a Groebner basis 
    with a set of variables. It is important that a proper elimination
    order is used to ensure that the variables that are to be eliminated
    by intersection can actually be eliminated.

    An example would be lexico-graphical ordering where the variables that are
    to be eliminated are the first ones. Or alternatively and superior in
    practice is to use block ordering where the variables that are to be
    eliminated are in the larger block w.r.t. to the ordering.

    # Arguments
    - `G`: a Groebner basis w.r.t a monomial odering that eliminates all 
       variables except the ones in `S_vars`
    - `S_vars`: the list of variables that are not to be eliminated

    # Returns
    - the intersection of the Groebner basis with the subring S
"""
function intersect(G, S_vars)
    sub_ideal = []
    for generator in G
        symbols = [Symbol(var) for var in AbstractAlgebra.vars(generator)]
        if issubset(symbols, S_vars)
            push!(sub_ideal, generator)
        end
    end
    return sub_ideal
end

"""
    _manage_rings(ideal, derivatives, R)

    This function creates a new ring that permits the elimination of variables.

    # Arguments
    - `ideal`: an ideal
    - `derivatives`: a dictionary of derivatives
    - `R`: a polynomial ring with two blocks of variables, where the first 
       block corresponds to the variables in `R` that are not in `derivatives` 
       and the second block corresponds to the variables in `derivatives`

    # Returns
    - the new ring
    - the new variables
    - a map from the old variables to the new variables
"""
function _manage_rings(derivatives, R)
    S_proper = []
    R_elim = []
    s = length(R.data.S) - length(derivatives)
    map = Dict()
    
    index_n_el = s + 1
    index_el = 1

    # Move the variables to the correct blocks and create a map
    for var in R.data.S
        if Symbol(var) in [Symbol(key) for key in keys(derivatives)]
            push!(S_proper, Symbol(var))
            map[Symbol(var)] = index_n_el
            index_n_el += 1
        else    
            push!(R_elim, Symbol(var))
            map[Symbol(var)] = index_el
            index_el += 1
        end
    end

    # Create new ring, where the first block is the eliminated variables
    R_vars_proper = append!(R_elim, S_proper)
    R_new_vars = [string(var) for var in R_vars_proper]
    R_new, R_new_vars = AlgebraicSolving.polynomial_ring(
        base_ring(R), R_new_vars, internal_ordering=:degrevlex)
    
    return R_new, R_new_vars, map
end


"""
    _swap_vars(poly, R_new, R_new_vars, map)

    This function swaps the variables in a polynomial to a new ring.

    # Arguments
    - `poly`: a polynomial
    - `R_new`: a new ring
    - `R_new_vars`: the new variables
    - `map`: a map from the old variables to the new variables

    # Returns
    - the polynomial in the new ring
"""
function _swap_vars(poly, R_vars, R_new_vars, map)
    vars_subst = [R_new_vars[map[Symbol(var)]] for var in R_vars]
    poly_new = poly(vars_subst...)
    return poly_new
end


"""
    differential_basis(ideal, derivatives, R, nf=false, info_level=0)

    This function computes the differential Gröbner basis of an ideal or 
    an approximation of a full differential ideal if not all derivatives
    are known. See Algorithm 4 in the thesis for reference.

    # Arguments
    - `ideal`: an ideal
    - `derivatives`: a dictionary of derivatives
    - `R`: a polynomial ring with two blocks of variables, where the first 
       block corresponds to the variables in `R` that are not in `derivatives` 
       and the second block corresponds to the variables in `derivatives`
    - `R_vars`: the variables in the ring, only needed for elimination
    - `nf`: a boolean indicating whether to compute the normal form
    - `info_level`: an integer indicating the level of information to print

    # Returns
    - the differential Gröbner basis
"""
function differential_basis(ideal, derivatives, R, R_vars = [], nf=false, info_level=0)
    n = R.data.nvars
    s = length(derivatives)    
    eliminate = n - s

    # If elimination is necessary reorganize the ring and substitute the variables
    if eliminate > 0
        R_new, R_new_vars, map = _manage_rings(derivatives, R)
        
        ideal_new_gens = [_swap_vars(elem, R_vars, R_new_vars, map) for elem in ideal.gens]
        ideal = AlgebraicSolving.Ideal(ideal_new_gens)

        derivatives_new = Dict()
        for (var, expr) in derivatives
            derivatives_new[_swap_vars(var, R_vars, R_new_vars, map)] = _swap_vars(expr, R_vars, R_new_vars, map)
        end
        derivatives = derivatives_new
    end 

    # Infer and create the subring
    S_vars = [Symbol(var) for var in R_new_vars[eliminate+1:end]]
    
    # Start computing the differential basis
    G1 = groebner_basis(ideal)
    pG1 = [diff_op(g, derivatives) for g in intersect(G1, S_vars)]
    if nf
        pG1 = [AlgebraicSolving.normal_form(pg, AlgebraicSolving.Ideal(G1)) for pg in pG1]
    end
    append!(pG1, G1)
    G2 = groebner_basis(AlgebraicSolving.Ideal(pG1), eliminate=eliminate,
                        intersect=false, info_level=info_level)
    if info_level > 0
        i = 1
        println("iteration ", i)
        println("#G = ", length(G1))
    end

    # Repeat until closed under diff_op
    while G1 != G2
        if info_level > 0
            i += 1
            println("iteration ", i)
            println("#G = ", length(G2))
        end
        G1 = G2
        pG1 = [diff_op(g, derivatives) for g in intersect(G1, S_vars)]
        if nf == true
            pG1 = [AlgebraicSolving.normal_form(pg, AlgebraicSolving.Ideal(G1)) for pg in pG1]
        end
            append!(pG1, G1)
        G2 = groebner_basis(AlgebraicSolving.Ideal(pG1), eliminate=eliminate,
                            intersect=false, info_level=info_level)
    end
    return G1
end