using AlgebraicSolving

"""
    chem_1()

    This function contains the implementation of a chemical system. 
    For reference, please see Example 1 in the following paper: 
    [link](https://arxiv.org/abs/1909.13608).

    The example can be found as example 2.5.4 in the thesis.
"""
function chem_1()
    
    R, (x1, x2, x3, x4, k1, k2, k3, T1, T2) = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),["x1","x2","x3","x4","k1","k2","k3","T1","T2"], 
        internal_ordering=:degrevlex)
    
    derivatives = Dict(
        x1 => - k1*x1 + k2*x2*x3,    
        x2 => k1*x1 - k2*x2*x3,
        x3 => - k2*x2*x3 + k3*x4,
        x4 => k2*x2*x3 - k3*x4,
        k1 => 0, k2 => 0, k3 => 0, T1 => 0, T2 => 0)

    ideal = AlgebraicSolving.Ideal([x1 + x2 - T1, x3 + x4 - T2])
    return ideal, derivatives, R
end

"""
    akzo_nobel()
 
    For reference, please see section 1.4 in the following paper: 
    [link](https://www.math.auckland.ac.nz/deptdb/dept_reports/497.pdf).

    The example can be found as example 2.5.5 in the thesis.
"""
function akzo_nobel()
    R, variables = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["dx6","dx7",
        "x1","x2","x3","x4","x5","x6","x7",
        "k1","k2","k3","k4","K","klA","Ks","pCO2","H"],
        internal_ordering=:degrevlex)

    # 1/H => H, 1/K => K
    (dx6, dx7, 
    x1, x2, x3, x4, x5, x6, x7, 
    k1, k2, k3, k4, K, klA, Ks, pCO2, H) = variables
    
    r1 = k1 * x1^4 * x7
    r2 = k2 * x3 * x5
    r3 = k2 * K * x1 * x5
    r4 = k3 * x1 * x4^2
    r5 = k4 * x6 * x7
    Fin = klA * (pCO2 * H - x2)

    derivatives = Dict(
        x1 => -2*r1 + r2 - r3 - r4,
        x2 => -r1 - 2*r4 - r5 + 2*Fin,
        x3 => -r1 - r2 + r3,
        x4 => -r2 + r3 - 2*r4,
        x5 => r2 - r3 + r5,
        x6 => dx6,
        x7 => dx7,
        k1 => 0,
        k2 => 0,
        k3 => 0,
        k4 => 0,
        K => 0,
        klA => 0,
        Ks => 0,
        pCO2 => 0,
        H => 0)

    ideal = AlgebraicSolving.Ideal([Ks * x1 * x4 - x6, x7^2 - x2])
    return ideal, derivatives, R
end

"""
    fast_slow_reaction()
 
    For reference, please see section 2.2 in the following paper: 
    [link](https://www.sciencedirect.com/science/article/pii/S1474667017433969).

    The example can be found as example 2.5.6 in the thesis.
"""
function fast_slow_reaction()
    R, variables = AlgebraicSolving.polynomial_ring(AlgebraicSolving.GF(101),
        ["dRA","dRB","dV2","dCAi","dF","dFi",
        "RA","RB","V1","V2","CA","CAi","CB","CC","F","Fi",
        "Keq","kB"],
        internal_ordering=:degrevlex)

    # 1/Keq => Keq
    (dRA, dRB, dV2, dCAi, dF, dFi, 
    RA, RB, V1, V2, CA, CAi, CB, CC, F, Fi, 
    Keq, kB) = variables
    
    derivatives = Dict(
        RA => dRA,
        RB => dRB,
        V1 => Fi - F,
        V2 => dV2,
        CA => Fi * V2 * (CAi - CA) - RA,
        CAi => dCAi,
        CB => - Fi * V2 * CB + RA - RB,
        CC => - Fi * V2 * CC + RB,
        F => dF,
        Fi => dFi,
        Keq => 0,
        kB => 0)

    ideal = AlgebraicSolving.Ideal([CA - CB * Keq, RB - kB * CB, V1*V2-1])
    return ideal, derivatives, R
end
