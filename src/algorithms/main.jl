using AbstractAlgebra
using AlgebraicSolving

"""
    _diff_op(q, derivatives)

    This function evaluates the linear differential operator from 
    Definition 7 and Notation 3 in the thesis for a polynomial and derivatives.
    The inputs should not explicitly depend on time.

    # Arguments
    - `q`: a polynomial
    - `derivatives`: a dictionary of derivatives

    # Returns
    - the linear differential operator applied to the polynomial
"""
function _diff_op(q, derivatives)    
    n = length(q.parent.data.S)
    @assert length(derivatives) <= n 
        "There can't be more derivatives than variables."
    
    result = 0
    for (var, value) in derivatives
        result += value * AbstractAlgebra.derivative(q, var)
    end
    return result
end

"""
    _intersect(G, S_vars)

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
function _intersect(G, S_vars)
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
    map_old_new = Dict()
    
    index_n_el = s + 1
    index_el = 1

    # Move the variables to the correct blocks and create a map
    for var in R.data.S
        if Symbol(var) in [Symbol(key) for key in keys(derivatives)]
            push!(S_proper, Symbol(var))
            map_old_new[Symbol(var)] = index_n_el
            index_n_el += 1
        else    
            push!(R_elim, Symbol(var))
            map_old_new[Symbol(var)] = index_el
            index_el += 1
        end
    end

    # Create new ring, where the first block is the eliminated variables
    R_vars_proper = append!(R_elim, S_proper)
    R_n_vars = [string(var) for var in R_vars_proper]
    R_n, R_n_vars = AlgebraicSolving.polynomial_ring(
        base_ring(R), R_n_vars, internal_ordering=:degrevlex)
    
    return R_n, R_n_vars, map_old_new
end


"""
    _swap_vars(poly, R_n, R_n_vars, map_old_new)

    This function swaps the variables in a polynomial to a new ring.

    # Arguments
    - `poly`: a polynomial
    - `R_n`: a new ring
    - `R_n_vars`: the new variables
    - `map_old_new`: a map from the old variables to the new variables

    # Returns
    - the polynomial in the new ring
"""
function _swap_vars(poly, R_vars, R_n_vars, map_old_new)
    vars_subst = [R_n_vars[map_old_new[Symbol(var)]] for var in R_vars]
    if poly == 0
        return poly
    end
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
function differential_basis(
        ideal, derivatives, R, R_vars = [], nf=false, info_level=0
    )
    
    n = R.data.nvars
    s = length(derivatives)    
    eliminate = n - s

    # If elimination is necessary reorganize the ring
    if eliminate > 0
        R_n, R_n_vars, map_old_new = _manage_rings(derivatives, R)
        
        ideal_new_gens = [
            _swap_vars(f, R_vars, R_n_vars, map_old_new) for f in ideal.gens]
        ideal = AlgebraicSolving.Ideal(ideal_new_gens)  

        derivatives_new = Dict()
        for (var, expr) in derivatives
            derivatives_new[_swap_vars(var, R_vars, R_n_vars, map_old_new)] = 
                _swap_vars(expr, R_vars, R_n_vars, map_old_new)
        end
        derivatives = derivatives_new

        # Infer and create the subring
        S_vars = [Symbol(var) for var in R_n_vars[eliminate+1:end]]
    else 
        S_vars = R.data.S
    end 

    # Start computing the differential basis
    G1 = groebner_basis(ideal)

    pG1 = [_diff_op(g, derivatives) for g in _intersect(G1, S_vars)]
    if nf
        reducer = AlgebraicSolving.Ideal(G1)
        pG1 = [AlgebraicSolving.normal_form(pg, reducer) for pg in pG1]
    end

    append!(pG1, G1)
    pG1 = Vector{typeof(ideal[1])}(pG1)
    G2 = groebner_basis(AlgebraicSolving.Ideal(pG1), eliminate=eliminate,
                        intersect=false, info_level=info_level)
    if info_level > 0
        i = 1
        println("iteration ", i)
        println("#G = ", length(G1))
    end

    # Repeat until closed under _diff_op
    while G1 != G2
        if info_level > 0
            i += 1
            println("iteration ", i)
            println("#G = ", length(G2))
        end
        G1 = G2
        pG1 = [_diff_op(g, derivatives) for g in _intersect(G1, S_vars)]
        if nf == true
            reducer = AlgebraicSolving.Ideal(G1)
            pG1 = [AlgebraicSolving.normal_form(pg, reducer) for pg in pG1]
        end
        
        append!(pG1, G1)
        pG1 = Vector{typeof(ideal[1])}(pG1)
        G2 = groebner_basis(AlgebraicSolving.Ideal(pG1), eliminate=eliminate,
                            intersect=false, info_level=info_level)
    end
    return G1
end