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
    R, variables = polynomial_ring(
        GF(101),["x6","x1","x2","x3","x4","x5","x6","x7","x8",
        "k1","k2","k3","k4","Kinv","klA","Ks","pCO2","Hinv"], 
        internal_ordering=:degrevlex)

        (x6, x1, x2, x3, x4, x5, x6, x7, x8, 
        k1, k2, k3, k4, Kinv, klA, Ks, pCO2, Hinv) = variables
    
    r1 = k1 * x1^4 * x7
    r2 = k2 * x3 * x5
    r3 = k2 * Kinv * x1 * x5
    r4 = k3 * x1 * x4^2
    r5 = k4 * x6 * x7
    Fin = klA * (pCO2 * Hinv - x2)

    derivatives = Dict(
        x1 => -2*r1 + r2 - r3 - r4,
        x2 => -r1 - 2*r4 - r5 + 2*Fin,
        x3 => -r1 - r2 + r3,
        x4 => -r2 + r3 - 2*r4,
        x5 => r2 - r3 + r5,
        x7 => x8,
        k1 => 0,
        k2 => 0,
        k3 => 0,
        k4 => 0,
        Kinv => 0,
        klA => 0,
        Ks => 0,
        pCO2 => 0,
        Hinv => 0)

    ideal = Ideal([Ks * x1 * x4 - x6, x7^2 - x2])
    return ideal, derivatives, R
end