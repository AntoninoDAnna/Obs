abstract type BC end

"""
        OBC
empty structure used to exploit julia automatic dispacth. It selects the function to use when Open Boundary Conditions applies
"""
struct OBC<:BC end
"""
        PBC
empty structure used to exploit julia automatic dispacth. It selects the function to use when Periodic Boundary Conditions applies
"""
struct PBC<:BC end

Base.length(::BC) = 1
Base.iterate(a::BC,state=1) = state>1 ? nothing : (a,state+1)
