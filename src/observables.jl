@doc raw"""
    mpcac(a0p,pp [,ca], bnd::Boundary = open)

Computes the effective PCAC mass associated to the correlators `a0p` and a `pp`
according to the boundary condition `bnd`. If `ca` is given, it first improve
the `a0p` correlator according to  `a0p-ca der_pp` and then computes the pcac mass.
It uses the symmetryc derivatives where necessary.

## Output specifics

If `a0p` and/or `pp` are `AbstractCorr`, then according to `bnd` the code uses:
   - `bnd == open`: the `obs` field where it strips the first and last component.
     This because in OBC  the first and last component lay on the boundary and are `0`
   - `bnd == periodic`: the `obs` field without modification

See also [`sym_der`](@ref), [`Boundary`](@ref), [`pa0_imp`](@ref)
"""
function mpcac(a0p,pp,bnd::B where {B<:BC})
    der_a0p = sym_der(a0p,bnd)
    return -der_a0p./(2 .*pp);
end

mpcac(a0p::T,pp::T,bnd::OBC) where T<:AbstractCorr = mpcac(a0p.obs[2:end-1], pp.obs[2:end-1],bnd)
mpcac(a0p::T,pp::T,bnd::PBC) where T<:AbstractCorr = mpcac(a0p.obs,pp.obs)
mpcac(a0p,pp,ca,bnd::B where {B<:BC} = OBC()) = mpcac(pa0_imp(a0p,pp,bnd,ca=ca),pp)
mpcac(a0p,pp) = mpcac(a0p,pp,OBC())

@doc """
        meff(v,x0,::B where {B<:BC})
        meff(v::T,::Union{OBC,OBC})

It computes the effective mass, using a consistent definition for both the forward and backward propagating correlator and assuming open boundary conditions applies
In particular, it uses the formula
```math
        m(t) = log(c(t)/c(t+1))
```
where `t=abs(src(v)-i)`.

the last parameter used to discriminate between Open and Periodic Boundary conditions.
"""
function meff(v,x0,::OBC)
    res = fill(one(eltype(v)),length(v)-2)
    for i in eachindex(res)
        res[i] = i< x0 ? log(abs(v[i+1]/v[i])) : log(abs(v[i]/v[i+1]))
    end
    return res
end

function meff(v,x0,::PBC)
    T = length(v)
    res = fill(one(eltype(v)),length(v))
    for i in eachindex(res)
        res[i] = i< x0 ? log(abs(v[i%T+1]/v[i])) : log(abs(v[i]/v[i%T+1]))
    end
    return res
end


meff(v,bnd::B where {B<:BC}) = meff(v,1,bnd)
meff(v::T,bnd::OBC) where T<:AbstractCorr = meff(v.obs[2:end-1], ObsIO.src(v),bnd)
meff(v::T,bnd::PBC) where T<:AbstractCorr = meff(v.obs, ObsIO.src(v),bnd)

meff(v) = meff(v,OBC())

@doc raw"""
    RI(HtoL::T, H::T, L::T, EL::T, EH::T) where T
    RI(HtoL::T, H::T, L::T, EL::T, EH::T) where T
    RI(HtoL::C, H::C, L::C, EL::T, EH::T) where {C<:AbstractCorr,T}

Compute the matrix element associated to the semileptonic decay `H -> L ℓ \nu` using the ratio
```math
    R_I(t_H,t_L) = \sqrt(2E_L) HtoL(t_H)/\sqrt(H(t_H) L(t_L)) \times  \exp(\frac{1}{2} (E_L t + E_H t_H))
```
It assume that `HtoL` is ordered according to `t_H`, i.e `xsrc` correspond to the meson H
"""
function RI(HtoL, H, L, EL, EH)
    ts = length(HtoL) -1
    res = Vector{eltype(HtoL)}(undef,length(HtoL))
    for tH in 0:ts
        tL = ts- tH
        res[tH+1] = sqrt(2EL) * abs(HtoL[tH+1]) * exp(0.5*EL*tL + 0.5*EH*tH ) / sqrt(abs(H[tH+1]*L[tL+1]))
    end
    return res
end

function RI(HtoL::Corr{3}, H::Corr{2}, L::Corr{2}, EL::R, EH::R) where R
    # We expect the correlators to have the following mass content:
    # HtoL : (h,l,l)
    # H: (h,l)
    # L: (l,l)
    l,_ = kappa(L)
    h, = let
        aux = kappa(H)
        aux[2] ==l || error("[RI] L and H do not have compatible masses")
        aux
    end
    isreverse  = let
        aux = kappa(HtoL)
        aux[2] == l || error("[RI] spectator quark in HtoL is not compatible with L")
        aux[2] == l || aux[3] == l || error("[RI] HtoL and L are not compatible")
        aux[1] == h || aux[3] == h || error("[RI] HtoL and H are not compatible")
        aux[1] == l
    end ## if true, we have reverse HtoL so that it reflect the H-meson time dependence
    T = length(H.obs)
    xsrc,xsnk = src(HtoL),snk(HtoL)
    isbackward = xsrc > xsnk
    isbackward && (xsrc,xsnk = xsnk,xsrc)
    _r = xsrc+1:xsnk+1 ## take into account the C -> julia index convention
    _HtoL = view(HtoL.obs,_r)
    _H = view(H.obs,_r)
    _L = view(L.obs,_r)
    ## from now on we are working with physical part of the correlators
    xor(isreverse, isbackward) && (_HtoL = reverse(_HtoL))
    !isbackward && return RI(_HtoL,_H,_L,EL,EH)
    return RI(_HtoL,reverse(_H),reverse(_L), EL, EH)
end

@doc raw"""
    RII(HtoL, H, L, EL, EH; xsink, xsrc)

Compute the matrix element associated to the semileptonic decay `H -> L` using the ratio
```math
    R_{II}^2(t) = 2E_L \frac{HtoL(t)HtoL(x_{sink}-t)}{H(x_{sink})L(x_{sink})}
```
"""
function RII(HtoL,H,L,EL,EH;xsnk::Int64,xsrc::Int64)
    R = HtoL[xsrc:xsnk].*HtoL[xsnk:-1:xsrc]
    R = (2*EL).*R./(H[xsnk]*L[xsnk]);
    return sqrt.(abs.(R))  ## if Btopi is V_0, then we have a minus sign due to Time reversal times Charge conjugation. Instead of bothering with the extra check, we take the absolute value y ya está
end

@doc raw"""
    RIII(HtoL, H, L, EL, EH; xsink, xsrc)

Compute the matrix element associated to the semileptonic decay `H -> L` using the ratio from 1903.05870 eq 5.1
```math
    R_{III}(t) = \sqrt(2E_L) HtoL(t)/\sqrt(H(x_{sink}) L(x_{sink})) \times  \exp(\frac{x_{sink}-2t}{2} (E_H - E_L))
```
"""
function RIII(HtoL,H, L, EL,EH;xsnk::Int64,xsrc::Int64)
    R = HtoL[xsrc:xsnk] ./sqrt(abs(H[xsnk]*L[xsnk]))
    aux2 = [exp(0.5*(EH-EL)*(xsnk-2*x+xsrc)) for x in xsrc:xsnk]
    return sqrt(2*EL).*abs.(R) .*aux2
end


@doc raw"""
     ps_dec(pa0,pp,mps,y0)

It computes the effective pseudoscalar decay constant associated to the correlator `pa0` and `pp`.
`mps` is the pseudoscalar mass and y0 is the source position.

If `pa0` and `pp` are a `AbstractCorr`, `y0` can be omitted and it will take the source position in `pa0`
"""
function ps_dec(pa0,pp,mps,y0)
    R= pa0./sqrt.(abs.(pp))
    e = exp.(0.5*abs.(mps).*[abs(i-y0) for i in 1:lastindex(pa0)])
    return sqrt.(2 ./abs.(mps)).*R.*e
end

function ps_dec(pa0::T,pp::T,mps,y0=0) where T<:AbstractCorr
    check_corr(pa0,pp,flag=no_gamma)
    if y0==0
        return ps_dec(pa0.obs[2:end-1],pp.obs[2:end-1],mps,pa0.y0)
    else
        return ps_dec(pa0.obs[2:end-1],pp.obs[2:end-1],mps,y0)
    end
end

@doc raw"""
     dec(c,m,y0)

It computes the effective decay constant associated to the correlator `c`.
`m` is the mass associated to it and `y0` the source position.

If `c` is a `AbstractCorr`, `y0` can be omitted and it will take the source position in `c`
"""
function dec(c,m,y0)
    e = exp.(0.5*abs.(m).*[abs(i-y0) for i in 1:lastindex(c)])
    return sqrt.(2 ./abs.(m)).*sqrt.(abs.(c)).*e
end

function dec(c::T,m,y0=0) where T<:AbstractCorr
    if y0==0
        return dec(c.obs[2:end-1],m,c.y0)
    else
        return dec(c.obs[2:end-1],m,c.y0)
    end
end
