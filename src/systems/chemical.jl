using AlgebraicSolving: polynomial_ring, GF, Ideal

function chem_1()
    # See for reference example 1 in 1909.13608
    R, (x1, x2, x3, x4, k1, k2, k3, T1, T2) = polynomial_ring(
        GF(101),["x1","x2","x3","x4","k1","k2","k3","T1","T2"], 
        internal_ordering=:degrevlex)
    derivatives = Dict(
        x1 => - k1*x1 + k2*x2*x3,    
        x2 => k1*x1 - k2*x2*x3,
        x3 => - k2*x2*x3 + k3*x4,
        x4 => k2*x2*x3 - k3*x4
    )
    ideal = Ideal([x1 + x2 - T1, x3 + x4 - T2])
    return ideal, derivatives, R
end

function akzo_nobel()

end