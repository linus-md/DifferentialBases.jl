using AlgebraicSolving: polynomial_ring, GF, Ideal

function simple_pendulum()
    R, (dl,l,v,u,y,x) = polynomial_ring(
        GF(101),["dl","l","v","u","y","x"], 
        internal_ordering=:degrevlex)
    derivatives = Dict(
        x => u,
        y => v,
        u => x*l,
        v => y*l - 1,
        l => dl)
    ideal = Ideal([x^2 + y^2 - 1])
    return ideal, derivatives, R
end

function masspoint_parabola()
    R, (dl, p1, p2, p3, v1, v2, v3, l) = polynomial_ring(
        GF(101),["dl","p1","p2","p3","v1","v2","v3","l"], 
        internal_ordering=:degrevlex)
    derivatives = Dict(
        p1 => v1,
        p2 => v2,
        p3 => v3,
        v1 => 2*l*p1,
        v2 => 2*l*p2,
        v3 => l - 1,
        l => dl)
    ideal = Ideal([p1^2 + p2^2 - p3])
    return ideal, derivatives, R
end