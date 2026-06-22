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
        m(t) = log(v(t)/v(t+1))
```
where `t=src(v)-i`.

the last parameter used to discriminate between Open and Periodic Boundary conditions.
"""
meff(v,::OBC) = log.(abs.(v[1:end-1]./v[2:end]))
meff(v::T,bnd::OBC) where T<:AbstractCorr = meff(v.obs[2:end-1],bnd)

meff(v,::PBC) = [log(abs(v[i]/v[i*length(v)+1])) for i in eachindex(v)]
meff(v::T,bnd::PBC) where T<:AbstractCorr = meff(v.obs,bnd)

meff(v,bnd::B where {B<:BC}) = meff(v,bnd)
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
    res = similar(HtoL)
    for tH in 0:ts
        tL = ts- tH
        res[tH+1] = sqrt(2EL) * abs(HtoL[tH+1]) * exp(0.5*EL*tL + 0.5*EH*tH ) / sqrt(abs(H[tH+1]*L[tL+1]))
    end
    return res
end

@doc raw"""
    RII(HtoL::T, H::T, L::T, EL::T, EH::T) where T
    RII(HtoL::T, H::T, L::T, EL::T, EH::T) where T
    RII(HtoL::C, H::C, L::C, EL::T, EH::T) where {C<:AbstractCorr,T}

Compute the matrix element associated to the semileptonic decay `H -> L ℓ \nu` using the ratio
```math
    R_II(t_H,t_L) = \sqrt(2E_L) \sqrt(HtoL(t_H) HtoL(t_L))/\sqrt(H(t_H + t_L) L(t_H +t_L))
```
It assume that `HtoL` is ordered according to `t_H`, i.e `xsrc` correspond to the meson H
"""
function RII(HtoL, H, L, EL, EH)
    R = HtoL.*reverse(HtoL)
    R = (2*EL).*R./(H[end]*L[end]);
    return sqrt.(abs.(R))
end

@doc raw"""
    RIII(HtoL::T, H::T, L::T, EL::T, EH::T) where T
    RIII(HtoL::T, H::T, L::T, EL::T, EH::T) where T
    RIII(HtoL::C, H::C, L::C, EL::T, EH::T) where {C<:AbstractCorr,T}

Compute the matrix element associated to the semileptonic decay `H -> L` using the ratio from 1903.05870 eq 5.1
```math
    R_{III}(t_H,t_L) = \sqrt(2E_L) HtoL(t_H)/\sqrt(H(t_H+t_L) L(t_H+t_L)) \times  \exp((t_H-t_L)/2 (E_H - E_L))
```
"""
function RIII(HtoL,H, L, EL,EH)
    R = HtoL ./sqrt(abs(H[end]*L[end]))
    ts = length(HtoL)-1
    aux2 = [exp(0.5*(EH-EL)*(2tH-ts)) for tH in 0:ts]
    return sqrt(2*EL).*abs.(R) .*aux2
end


@doc raw"
     __deduce_mass(HtoL, H, L)

read the quark masses in the correlators and deduce the quark contents.
It returns a tuple with `(h,sp,l)` where `h` is the heavy quark mass, `sp` is the specator quark and `l` is the light quark.
"
function __deduce_mass(HtoL::Corr{3}, H::Corr{2}, L::Corr{2})
    function get_sp(kappa, (h1,h2), (l1,l2))
        for k in kappa
            if ((k == h1) || (k == h2)) && ((k == l1) || (k==l2))
               return  k
            end
        end
        error("Cannot deduce spectator quark")
    end
    function get_other((k1,k2),sp)
        k1 == sp && return k2
        k2 == sp && return k1
        error("correlator with masses $k1, $k2 does not have the spectators quark $sp")
    end
    k3,kh,kl =  kappa(HtoL), kappa(H), kappa(L)
    sp = get_sp(k3,kh,kl)
    h = get_other(kh,sp)
    l = get_other(kl,sp)
    return h,sp,l
end

function _ratios(F::Function, HtoL::Corr{3}, H::Corr{2},
                L::Corr{2}, EL::R, EH::R) where R
    # We expect the correlators to have the following mass content:
    # HtoL : (h,ls,l) where ls is the spectator quark
    # H: (h,ls)
    # L: (ls,l)
    h,sp,l = __deduce_mass(HtoL,H,L)
    isreverse  = let
        _h,_sp,_l = kappa(HtoL) ## we expect the 3pt to have these quarks
        _sp == sp || error("[$F] The spectator quark in HtoL is not in the correct position")
        _h == l || _l == l || error("[$F] HtoL and L are not compatible")
        _h == h || _l == h || error("[$F] HtoL and H are not compatible")
        _h == l
    end ## if true, we have reverse HtoL so that it reflect the H-meson time dependence
    T = length(H.obs)
    xsrc,xsnk = src(HtoL),snk(HtoL)
    isbackward = xsrc > xsnk
    if isbackward
        xsrc,xsnk = xsnk,xsrc
    end
    _r = xsrc+1:xsnk+1 ## take into account the C -> julia index convention
    _HtoL = HtoL.obs[_r]
    _H = H.obs[_r]
    _L = L.obs[_r]
    xor(isreverse, isbackward) && ( _HtoL = reverse(_HtoL))
    !isbackward && return F(_HtoL,_H,_L,EL,EH)
    return F(_HtoL,reverse(_H),reverse(_L), EL, EH)
end

RI(HtoL::Corr{3}, H::Corr{2}, L::Corr{2}, EL::R, EH::R) where R =
    _ratios(RI,HtoL,H,L,EL,EH)
RII(HtoL::Corr{3}, H::Corr{2}, L::Corr{2}, EL::R, EH::R) where R =
    _ratios(RII,HtoL,H,L,EL,EH)
RIII(HtoL::Corr{3}, H::Corr{2}, L::Corr{2}, EL::R, EH::R) where R =
    _ratios(RIII,HtoL,H,L,EL,EH)


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
