@doc raw"""
     v_imp(vv,vt...;cv, L=1, theta1,theta2, bnd::Boundary=open)
     v_imp(vv::AbstractCorr,vt::AbstractCorr...;cv,L=1, theta1=Float64[],theta2=Float64[], bnd::Boundary=open)

Compute the improved G_{VV} =  \sum_{k=1}^3G_{V_iV_i} correlator according to the
equations
```math
    G^I_{V_iV_i} = G_{V_i V_i} + 2ac_V \de_t G_{V_i T_i0} - 2c_V \sum_{j=1}^{3} i \sin(p^j a) G_{V_i T_ij}
```

# Arguments
  - vv: correlator G_{VV}. It has to be already averaged over the three polarizations.
  - vt: correlators G_{VT}. It should be either G_{ViT0i} (already averaged
        over the three polarizations) if the correlator is at `0` momentum or
        G_{ViT0i}, G_{V1T12},G_{V1T13}, G_{V2T12},G_{V2T23}, G_{V3T13},G_{V3T23} if at
        non-zero momentum

# Keyword Arguments
  - `cv`: mandatory keyword argument. It is the improvement coefficient
  - `L::Int64=1`: Lattice size. Needed to convert the theta angle into momenta.
  - `theta1=Float64[]`: theta angle of the first quark. If empty, is read from `vv`
  - `theta2=Float64[]`: theta angles of the second quark. If empty, is read from `vv`
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function v_imp(vivi,vit0i,bnd::OBC; cv)
    T = eltype(vivi)
    der_vt = [zero(T);sym_der(vit0i[2:end-1],bnd);zero(T)]
    return vivi.-2cv*der_vt
end

function v_imp(vivi,vit0i,bnd::B where {B<:BC}; cv)
    der_vt = sym_der(vit0i,bnd)
    return vivi.-2cv*der_vt
end

function v_imp(vivi, vit0i, v1t12, v1t13, v2t12, v2t23, v3t13, v3t23,bnd::B where {B<:BC};
               cv,L=1,theta1,thet2)
    p = (theta1.-theta2)/L
    imp = v_imp(vivi,vit0i,bnd,cv=cv)
    all(pp==0 for pp in p) && return imp
    aux = -sin(p[1]).*(v2t12.+v3t13).+
        sin(p[2]).*(v1t12.-v3t23).+
        sin(p[3]).*(v1t13.+v2t23)
    return imp .-2cv*aux
end

function v_imp(vivi::T, vit0i::T, bnd::OBC;cv)::T where T<: AbstractCorr
    imp = v_imp(vivi.obs[2:end-1],vit0i.obs[2:end-1],bnd,cv=cv)
    return ObsIO.__update__(vv,obs=imp)
end

function v_imp(vivi::T, vit0i::T, bnd::PBC;cv)::T where T<: AbstractCorr
    imp =v_imp(vivi.obs,vit0i.obs,bnd,cv=cv)
    return ObsIO.__update__(vv,obs=imp)
end

function v_imp(vivi::T, vit0i::T, v1t12::T,
               v1t13::T, v2t12::T, v2t23::T,
               v3t13::T, v3t23::T,bnd::OBC;
               cv,L=1,theta1=nothing,theta2=nothing)::T where T<: AbstractCorr

    isnothing(theta1) && (theta1 = theta(vv)[1])
    isnothing(theta2) && (theta2 = theta(vv)[2])
    imp = v_imp(vivi.obs, vit0i.obs,
                v1t12.obs,v1t13.obs,
                v2t12.obs,v2t23.obs,
                v3t13.obs,v3t23.obs,bnd,
                cv=cv,L=L,theta1=theta1,theta2=theta2)
    return ObsIO.__update__(vv,obs=imp)
end

function v_imp(vivi::T, vit0i::T, v1t12::T,
               v1t13::T, v2t12::T, v2t23::T,
               v3t13::T, v3t23::T,bnd::PBC;
               cv,L=1,theta1=nothing,theta2=nothing)::T where T<: AbstractCorr

    isnothing(theta1) && (theta1 = theta(vv)[1])
    isnothing(theta2) && (theta2 = theta(vv)[2])
    imp = v_imp(vivi.obs, vit0i.obs,
                v1t12.obs,v1t13.obs,
                v2t12.obs,v2t23.obs,
                v3t13.obs,v3t23.obs,bnd,
                cv=cv,L=L,theta1=theta1,theta2=theta2)
    return ObsIO.__update__(vv,obs=imp)
end

@doc raw"""
    pa0_imp(pa0,pp; ca)
    pa0_imp(pa0::T,pp::T, bnd::Union{OBC, PBC}; ca)::T where T<:AbstractCorr

Improve the correlator G_{PA_0} according to the equation
```math
    G_{PA_0^{I}} = G_{PA_0} - a c_A \\de_t G_{PP}
```

`bnd` is an expected parameter used to discriminate between Open and Periodic Boundary Contidions
# Keyword Arguments
  - `ca`: mandatory keyword argument. It is the improved operator [See also `ca_fit`]

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function pa0_imp(pa0, pp,bnd::OBC; ca)
    T = eltype(pa0)
    der_p= [zero(T);sym_det(pp[2:end-1],bnd);zero(T)]
    return pa0.-ca.*der_p
end

function pa0_imp(pa0, pp,bnd::B where {B<:BC}; ca)
    der_p=sym_der(pp,bnd)
    return pa0.-ca.*der_p
end


function pa0_imp(pa0::T, pp::T, bnd::OBC; ca)::T where T<:AbstractCorr
    imp = pa0_imp(pa0.obs,pp.obs,bnd,ca=ca)
    return ObsIO.__update__(pa0,obs=[_zero;imp;_zero])
end

function pa0_imp(pa0::T,pp::T,bnd::PBC;ca)::T where T<:AbstractCorr
    imp = pa0_imp(pa0.obs,pp.obs,bnd,ca=ca,)
    return ObsIO.__update__(pa0,obs=imp)
end

@doc """
    a0a0_imp(a0a0,pa0; ca,  bnd::Boundary=open)
    a0a0_imp(a0a0::AbstractCorr,pa0::AbstractCorr; ca,  bnd::Boundary=open)::AbstractCorr

Improve the correlator G_{A_0 A_0} according to the equation
```math
    G_{A_0^{I} A_0^{I}} = G_{A_0 A_0} - 2 a c_A \\de_t G_{P A_0}
```
# Keyword Arguments
  - `ca`: mandatory keyword argument. It is the improved operator [See also `ca_fit`]
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function a0a0_imp(a0a0, pa0, bnd::OBC; ca)
    T = eltype(a0a0)
    der_pa0 = [zero(T);sym_der(pa0[2:end-1],bnd); zero(T)]
    return a0a0 .- (2*ca).*der_pa0
end

function a0a0_imp(a0a0, pa0, bnd::B where {B<:BC}; ca)
    der_pa0 = sym_der(pa0,bnd)
    return a0a0 .- (2*ca).*der_pa0
end

function a0a0_imp(a0a0::T, pa0::T,bnd::OBC; ca )::T where T<:AbstractCorr
    imp = a0a0_imp(a0a0.obs,pa0.obs,bnd,ca=ca)
    return ObsIO.__update__(a0a0,obs=imp)
end

function a0a0_imp(a0a0::T, pa0::T,bnd::PBC; ca )::T where T<:AbstractCorr
    imp = a0a0_imp(a0a0.obs,pa0.obs,bnd,ca=ca)
    return ObsIO.__update__(a0a0,obs=imp)
end

@doc """
     pv_imp(pv,pt...;cv,L::Int64=1,theta1, theta2, bnd::Boundary = open)
     pv_imp(pv::T, pt::T...; cv, L=1, theta1=Float64[],theta2 = Float64[], bnd::Boundary=open)

Compute the improved G_{PV} =  1/3 \\sum_{k=1}^3G_{PV_i} correlator according to the equations
```math
  G^I_{PV_i} = G_{P V_i} + ac_V \\de_t G_{P T_i0} - c_V \\sum_{j=1}^{3} \\sin(p^j a) G_{P T_ij}
```

# Arguments
  - pv: correlator G_{PV}. It has to be already averaged over the three polarizations.
  - pt: correlators G_{PT}. It should be either G_{PT0i} (already averaged
        over the three polarizations) if the correlator is at `0` momentum or
        G_{PT0i}, G_{PT12},G_{PT13},G_{PT23} if at non-zero momentum

# Keyword Arguments
  - `cv`: mandatory keyword argument. It is the improvement coefficient
  - `L::Int64=1`: Lattice size. Needed to convert the theta angle into momenta.
  - `theta1=Float64[]`: theta angle of the first quark. If empty, is read from `pv`
  - `theta2=Float64[]`: theta angles of the second quark. If empty, is read from `pv`
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives
  - `real`: mandatory flag that indicate whether we are improving the real or imaginary part of ´G_{PVi}´. WARNING to improve the real (imaginary) part of ´G_{PVi}´ it is the real (imaginary) part of ´G_{PT0i}´  and the imaginary (real) part of ´G_{PTij` are required.
See also [`sym_der`](@ref), [`Boundary`](@ref)
            """

function pv_imp(pvi, pt0i,bnd::OBC; cv,real::Bool)
    T = eltype(pvi)
    der_pt = [zero(T);sym_der(pt0i[2:end-1],bnd); zero(T)]
    return pvi .- cv.*der_pt # we only have access to T_{0i}, but we want T_{i0}
end

function pv_imp(pvi, pt0i,bnd::B where {B<:BC};cv,real::Bool)
    der_pt = sym_der(pt0i,bnd)
    return pvi .- cv.*der_pt # we only have access to T_{0i}, but we want T_{i0}
end

function pv_imp(pv,pt0i,pt12,pt13,pt23,bnd::B where {B<:BC}; cv,L::Int64=1,theta1, theta2,real::Bool)
    imp = pv_imp(pv.pt01,bnd,cv=cv,real=real)
    p = (theta1.-theta2)./L
    (all(p.==0) || all(p[2:end].==p[1]) ) && return imp
    aux  = -sin(p[1])*(pt12+pt13)
    aux +=  sin(p[2])*(pt12-pt23)
    aux +=  sin(p[3])*(pt13+pt23)
    return real ? imp .+ cv.*aux./3 : imp .- cv.*aux./3
end

function pv_imp(pvi::T, pt0i::T, bnd::OBC; cv,real::Bool)::T where T<:AbstractCorr
    imp = pv_imp(pvi.obs,pt0i.obs,bnd,cv=cv,real=real)
    return ObsIO.__update__(pvi,obs=imp)
end

function pv_imp(pvi::T, pt0i::T, bnd::PBC; cv,real::Bool)::T where T<:AbstractCorr
    imp = pv_imp(pvi.obs,pt0i.obs,bnd,cv=cv,real=real)
    return ObsIO.__update__(pv,obs=imp)
end

function pv_imp(pv::T, pt0i::T, pt12::T, pt13::T, pt23::T, bnd::OBC; cv, L=1, theta1=nothing,theta2=nothing,real::Bool)::T where T<:AbstractCorr
    isnothing(theta1) && (theta1 = theta(pv)[1])
    isnothing(theta2) && (theta2 = theta(pv)[2])
    imp =  pv_imp(pv.obs, pt0i.obs, pt12.obs,
                  pt13.obs, pt23.obs,bnd,
                  cv=cv,L=L, theta1=theta1, theta2=theta2,real=real)
    return ObsIO.__update__(pv,obs=imp)
end

function pv_imp(pv::T, pt0i::T, pt12::T, pt13::T, pt23::T, bnd::PBC; cv, L=1, theta1=nothing,theta2=nothing, real::Bool)::T where T<:AbstractCorr
    isnothing(theta1) && (theta1 = theta(pv)[1])
    isnothing(theta2) && (theta2 = theta(pv)[2])
    imp =  pv_imp(pv.obs, pt0i.obs, pt12.obs, pt13.obs, pt23.obs,bnd,
                  cv=cv,L=L, theta1=theta1, theta2=theta2,real=real)
    return ObsIO.__update__(pv,obs=imp)
end

@doc raw"""
     pv0_imp(pv0, pt...; theta1,theta2, cv,L::Int64=1, bnd::Boundary=open,real::Bool)
     pv0_imp(pv0::T,pt::T ...;cv,L=1,theta1 = Float64[], theta2 = Float64[], bnd::Boundary=open,real::Bool)::T

Improve the correlator ´G_{PV0}´ and ´G_{PV0P}´ with the tensor correlator according to the improvement equation

```math
    V_0^I(t,\vec p) = V_0(t,\vec p) - c_V \sum_{k=1}^3 i sin(ap^k) T_{0k}(t,\vec p)
```

# Arguments
  - pv0: correlator ´G_{PV0}´
  - pt: correlators ´G_{PT}´. It contains the correlator ´G_{PT01}´, ´G_{PT02}´, ´G_{PT03}´

# Keyword Arguments
  - `cv`: mandatory keyword argument. It is the improvement coefficient
  - `L::Int64=1`: Lattice size. Needed to convert the theta angle into momenta.
  - `theta1=Float64[]`: theta angle of the first quark. If empty, is read from `pv`
  - `theta2=Float64[]`: theta angles of the second quark. If empty, is read from `pv`
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives
  - `real::Bool` : mandatory flag that indicate whether we are improving the real or imaginary part of ´G_{PV0P}´. WARNING to improve the real (imaginary) part of ´G_{PV0P}´ the imaginary (real) part of ´G_{PT0iP}´ is required

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function pv0_imp(pv0, pt01, pt02, pt03; theta1,theta2, cv,L::Int64=1,real::Bool)
    p = (theta1.-theta2)./L
    all(p.==0) &&   return pv0
    aux = sin(p[1]).*pt01.+sin(p[2]).*pt02.+sin(p[3]).*pt03

    return real ? pv0.+cv.*aux : pv0.-cv.*aux
end

function pv0_imp(pv0::T, pt01::T, pt02::T, pt03::T, ::OBC;cv,L=1,theta1 = nothing, theta2 = nothing,real::Bool)::T where T<:AbstractCorr
    isnothing(theta1) && (theta1 = theta(pv0)[1])
    isnothing(theta2) && (theta2 = theta(pv0)[2])
    imp =  pv0_imp(pv0.obs, pt01.obs, pt02.obs,
                   pt03.obs, cv=cv, L=L,theta1=theta1, theta2=theta2,real = real)
    return ObsIO.__update__(pv0,obs=imp)
end


function pv0_imp(pv0::T, pt01::T, pt02::T, pt03::T, ::PBC;cv,L=1,theta1 = nothing, theta2 = nothing,real::Bool)::T where T<:AbstractCorr
    isnothing(theta1) && (theta1 = theta(pv)[1])
    isnothing(theta2) && (theta2 = theta(pv)[2])
    imp =  pv0_imp(pv0.obs, pt01.obs, pt02.obs,
                   pt03.obs, cv=cv, L=L,theta1=theta1, theta2=theta2,real=real)
    return ObsIO.__update__(pv0,obs=imp)
end
