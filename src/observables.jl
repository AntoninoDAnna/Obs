@doc raw"""
    mpcac(a0p,pp [,ca], bnd::Boundary = open)

Computes the effective PCAC mass associated to the correlators `a0p` and a `pp`
according to the boundary condition `bnd`. If `ca` is given, it first improve
the `a0p` correlator according to  `a0p-ca der_pp` and then computes the pcac mass.
It uses the symmetryc derivatives where necessary.

## Output specifics

If `a0p` and/or `pp` are `juobs.Corr`, then according to `bnd` the code uses:
   - `bnd == open`: the `obs` field where it strips the first and last component.
     This because in OBC  the first and last component lay on the boundary and are `0`
   - `bnd == periodic`: the `obs` field without modification

See also [`sym_der`](@ref), [`Boundary`](@ref), [`pa0_imp`](@ref)
"""
function mpcac(a0p,pp,bnd::Boundary=open)
    der_a0p = sym_der(a0p,bnd)
    return -der_a0p./(2 .*pp);
end

function mpcac(a0p::juobs.Corr,pp::juobs.Corr,bnd::Boundary = open)
    if bnd == open
        return mpcac(a0p.obs[2:end-1], pp.obs[2:end-1],bnd)
    elseif bnd == periodic
        return mpcac(a0p.obs,pp.obs,bnd)
    end
end

mpcac(a0p,pp,ca,bnd::Boundary=open) = mpcac(pa0_imp(a0p,pp,ca=ca, bnd=bnd),pp, bnd)

@doc raw"""
     meff(v,bnd::Boundary=open)

Compute the effective mass associated to the correlator `v` according to the boundary condition
`bnd`.

##  Output Specifics

The length of the retuned vector is `length(v)-1` if OBC, while is `length(v)` if PBC

"""
function meff(v,bnd::Boundary=open)
    if bnd == open
        return log.(abs.(v[1:end-1]./v[2:end]));
    elseif bnd == periodic
        T = length(v)
        return log.([v[t]/v[t%T+1] for t in 1:T])
    end
end

function meff(v::juobs.Corr,bnd::Boundary=open)
    if bnd == open
        return meff(v.obs[2:end-1],bnd)
    elseif bnd == periodic
        return meff(v.obs,bnd)
    end
end

@doc raw"""
    RI(HtoL, H, L, EL, EH; xsink, xsrc)

Compute the matrix element associated to the semileptonic decay `H -> L` using the ratio
```math
    R_I(t) = \sqrt(2E_L) HtoL(t)/\sqrt(H(x_{sink}-t) L(t)) \times  \exp(\frac{1}{2} (E_L t + E_H(x_{sink}-t)))
```
"""
function RI(HtoL,H,L,EL,EH;xsnk::Int64,xsrc::Int64)
    aux = H[xsnk:-1:xsrc].*L[xsrc:xsnk]
    R   = abs.(HtoL[xsrc:xsnk])./sqrt.(abs.(aux));
    return sqrt(2EL).*R .*exp.(0.5.*[EL*(z0-xsrc)+EH*(xsnk-z0) for z0 in xsrc:xsnk])
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
    R = HtoL[xsrc:xsnk] ./sqrt(H[xsnk]*L[xsnk])
    aux2 = [exp(0.5*(EH-EL)*(xsnk-2*x+xsrc)) for x in xsrc:xsnk]
    return sqrt(2*EL).*abs.(R) .*aux2
end

for f in (:RI, :RII, :RIII)
    @eval function $(f)(HtoL::juobs.Corr,H::juobs.Corr,L::juobs.Corr,EL,EH;xsnk::Int64)
        if HtoL.theta1 != L.theta1
            error("[$f] theta mismatch. thetas: $(HtoL.theta1) and $(L.theta1)")
        end
        if HtoL.kappa[2] != H.kappa[1]
            error("[$f] heavy mass mismatch.  masses $(HtoL.kappa[2]) and $(H.kappa[1])")
        end
        if HtoL.kappa[1] != H.kappa[2] !=L.kappa[1]
            error("[$f] light mass mismatch.  masses $(HtoL.kappa[2]), $(H.kappa[1]) and $(L.kappa[1])")
        end
        if HtoL.y0 != H.y0 != L.y0
            error("[$f] source position mismatch. sources at $(HtoL.y0), $(H.y0) and $(L.y0)")
        end
        return $(f)(HtoL.obs[2:end-1],H.obs[2:end-1],L.obs[2:end-1],EL,EH,xsnk = xsnk, xsrc = HtoL.y0)
    end
end

@doc raw"""
     ps_dec(pa0,pp,mps,y0)

It computes the effective pseudoscalar decay constant associated to the correlator `pa0` and `pp`.
`mps` is the pseudoscalar mass and y0 is the source position.

If `pa0` and `pp` are a `juobs.Corr`, `y0` can be omitted and it will take the source position in `pa0`
"""
function ps_dec(pa0,pp,mps,y0)
    R= pa0./sqrt.(abs.(pp))
    e = exp.(0.5*mps.*[i-y0 for i in 1:lastindex(pa0)])
    return (sqrt(2/mps)).*R.*e
end

function ps_dec(pa0::juobs.Corr,pp::juobs.Corr,mps,y0=0)
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

If `c` is a `juobs.Corr`, `y0` can be omitted and it will take the source position in `c`
"""
function dec(c,m,y0)
    e = exp.(0.5*m.*[i-y0 for i in 1:lastindex(c)])
    return sqrt(2/m).*sqrt.(abs.(c)).*e
end

function dec(c::juobs.Corr,m,y0=0)
    if y0==0
        return dec(vv.obs[2:end-1],mvec,vv.y0)
    else
        return dec(vv.obs[2:end-1],mvec,y0)
    end
end
