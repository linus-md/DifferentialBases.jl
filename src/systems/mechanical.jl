using AlgebraicSolving: polynomial_ring, GF, Ideal
"""
    simple_pendulum()

    This function contains the implementation of a simple pendulum.
    For reference see Example 2.5.1 in the thesis.
"""
function simple_pendulum()
    R, variables = polynomial_ring(
        GF(101),["dl","x","y","u","v",""], 
        internal_ordering=:degrevlex)
    (dl,x,y,u,v,l) = variables

    derivatives = Dict(
        x => u,
        y => v,
        u => x*l,
        v => y*l - 1,
        l => dl)

    ideal = Ideal([x^2 + y^2 - 1])
    return ideal, derivatives, R
end

"""
    double_pendulum()

    This function contains the implementation of a double pendulum.
    For reference see Example 2.5.2 in the thesis.
"""
function double_pendulum()
    R, variables = polynomial_ring(GF(101),
        ["dl1","dl2","x1","y1","u1","v1","x2","y2","u2","v2","l1","l2"],
        internal_ordering=:degrevlex)
    (dl1,dl2,x1,y1,u1,v1,x2,y2,u2,v2,l1,l2) = variables

    derivatives = Dict(
        x1 => u1,
        y1 => v1,
        u1 => - l1*x1 - l2*(x1 - x2),
        v1 => - l1*y1 - l2*(y1 - y2) - 1,
        x2 => u2,
        y2 => v2,
        u2 => - l2*(x2 - x1),
        v2 => - l2*(y2 - y1) - 1,
        l1 => dl1,
        l2 => dl2)

    ideal = Ideal([x1^2 + y1^2 - 1, (x2-x1)^2 + (y2-y1)^2 - 1])
    return ideal, derivatives, R
end

"""
    triple_pendulum()

    This function contains the implementation of a triple pendulum.
    
    For reference see Example 2.5.2 in the thesis.
"""
function triple_pendulum()
    R, variables = polynomial_ring(
        GF(101),["dl1","dl2","dl3","x1","y1","u1","v1","x2","y2","u2","v2",
        "x3","y3","u3","v3","l1","l2","l3"],
        internal_ordering=:degrevlex)
    (dl1,dl2,dl3,x1,y1,u1,v1,x2,y2,u2,v2,x3,y3,u3,v3,l1,l2,l3) = variables
    
    derivatives = Dict(
        x1 => u1,
        y1 => v1,
        u1 => - l1*x1 - l2*(x1 - x2),
        v1 => - l1*y1 - l2*(y1 - y2) - 1,
        x2 => u2,
        y2 => v2,
        u2 => - l2*(x2 - x1) - l3*(x2 - x3),
        v2 => - l2*(y2 - y1) - l3*(y2 - y3) - 1,
        x3 => u3,
        y3 => v3,
        u3 => - l3*(x3 - x2),
        v3 => - l3*(y3 - y2) - 1,
        l1 => dl1,
        l2 => dl2,
        l3 => dl3)

    ideal = Ideal([x1^2 + y1^2 - 1, (x2-x1)^2 + (y2-y1)^2 - 1,
            (x3-x2)^2 + (y3-y2)^2 - 1])
    return ideal, derivatives, R
end

"""
    masspoint_parabola()

    This function contains the implementation of a masspoint restricted to
    a parabola.
    
    For reference see Example 2.5.3 in the thesis.
"""
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
