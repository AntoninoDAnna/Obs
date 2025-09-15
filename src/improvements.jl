import juobs, ADerrors

@doc raw"""
     v_imp(vv,vt...;cv, L=1, theta1,theta2, bnd::Boundary=open)
     v_imp(vv::juobs.Corr,vt::juobs.Corr...;cv,L=1, theta1=Float64[],theta2=Float64[], bnd::Boundary=open)

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
function v_imp(vv,vt...;cv, L=1, theta1,theta2, bnd::Boundary=open)
    der_vt = sym_der(vt[1],bnd)
    imp = vv.-2cv*der_vt

    p = (theta1.-theta2)/L
    if all(p.==0)
        return imp
    end

    if length(vt) !=7
        error("[v_imp]: unexpected number of tensor current correlators")
    end

    v1t12, v1t13, v2t12,v2t23, v3t13,v3t23 = [v for vt in v[2:end]]
    aux = -sin(p[1]).*(v2t12.+v3t13).+sin(p[2]).*(v1t12.-v3t23)+sin(p[3]).*(v1t13+v2t23)

    return imp .-2cv*aux
end

function v_imp(vv::juobs.Corr,vt::juobs.Corr...;cv,L=1, theta1=Float64[],theta2=Float64[], bnd::Boundary=open)
    check_corr(vv,vt..., flag = no_gamma)
    vv, vt = if bnd == open
        vv.obs[2:end-1], [x.obs[2:end-1] for x in vt]
        elseif bnd ==period
        vv.obs , getfield.(vt,:obs);
    end
    if isempty(theta1)
        theta1 = vv.theta1
    end
    if isempty(theta2)
        theta2 = vv.theta2
    end
    imp = v_imp(vv,vt..., cv = cv, L=L, theta=theta1, theta2=theta2)
    return juobs.Corr(imp, vv.kappa,vv.mu, vv.gamma,vv.y0,vv.theta1,vv.theta2)
end

@doc raw"""
    pa0_imp(pa0,pp; ca, bnd::Boundary =open)
    pa0_imp(pa0::juobs.Corr,pp::juobs.Corr; ca, bnd::Boundary =open)::juobs.Corr

Improve the correlator G_{PA_0} according to the equation
```math
    G_{PA_0^{I}} = G_{PA_0} - a c_A \\de_t G_{PP}
```
# Keyword Arguments
  - `ca`: mandatory keyword argument. It is the improved operator [See also `ca_fit`]
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function pa0_imp(pa0, pp; ca, bdn::Boundary=open)
    der_p=sym_der(pp,bnd)
    return corr_pa0.-ca.*der_p
end

function pa0_imp(pa0::juobs.Corr, pp::juobs.Corr; ca,  bdn::Boundary=open)
    check_corr(pa0, pp,flag=no_gamma)
    imp = if bnd == open
        pa0_imp(pa0.obs[2:end-1],pp.obs[2:end-1],ca=ca, bnd=bnd)
    elseif bnd==periodic
        pa0_imp(pa0.obs,pp.obs,ca=ca, bnd=bnd)
    end
    return juobs.Corr(imp,pp.kappa,pp.mu,pp.gamma, pp.y0,pp.theta1,pp.theta2)
end

@doc """
    a0a0_imp(a0a0,pa0; ca,  bdn::Boundary=open)
    a0a0_imp(a0a0::juobs.Corr,pa0::juobs.Corr; ca,  bdn::Boundary=open)::juobs.Corr

Improve the correlator G_{A_0 A_0} according to the equation
```math
    G_{A_0^{I} A_0^{I}} = G_{A_0 A_0} - 2 a c_A \\de_t G_{P A_0}
```
# Keyword Arguments
  - `ca`: mandatory keyword argument. It is the improved operator [See also `ca_fit`]
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function a0a0_imp(a0a0, pa0; ca, bdn::Boundary=open)
    der_pa0 = sym_der(pa0,bnd)
    return a0a0 .- (2*ca).*der_pa0
end

function a0a0_imp(a0a0::juobs.Corr, pa0::juobs.Corr; ca,  bdn::Boundary=open)
    check_corr(a0a0, pa0,no_gamma)
    imp = if bnd == open
        a0a0_imp(a0a0.obs[2:end-1],pa0.obs[2:end-1],ca=ca, bnd=bnd)
    elseif bnd==periodic
        a0a0_imp(a0a0.obs,pa0.obs,ca=ca, bnd=bnd)
        return juobs.Corr(imp,pa0.kappa,pa0.mu,pa0.gamma, pa0.y0,pa0.theta1,pa0.theta2)
    end
end

@doc """
     pv_imp(pv,pt...;cv,L::Int64=1,theta1, theta2, bnd::Boundary = open)
     pv_imp(pv::juobs.Corr, pt::juobs.Corr...; cv, L=1, theta1=Float64[],theta2 = Float64[], bnd::Boundary=open)

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

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function pv_imp(pv,pt...;cv,L::Int64=1,theta1, theta2, bnd::Boundary = open)
    der_pt = der_sym(pt[1])

    p = (theta1.-theta2)./L
    if all(p.==0) || all(p[2:end].==p[1])
        return pv .- cv.*der_pt # we only have access to T_{0i}, but we want T_{i0}
    end
    if length(pt) .!=4
        error("[pv_imp]: unexpected number of tensor current correlators")
    end

    pt12, pt13, pt23 = pt[2:end]
    aux  = -sin(p[1])*(pt12+pt13)
    aux +=  sin(p[2])*(p12-p23)
    aux +=  sin(p[3])*(pt13-pt23)

    return pv .-cv *der_vt .- cv*aux./3
end

function pv_imp(pv::juobs.Corr, pt::juobs.Corr...; cv, L=1, theta1=Float64[],theta2 = Float64[], bnd::Boundary=open)
    check_corr(pv,pt...,no_gamma)
    if isempty(theta1)
        theta1 = pv.theta1
    end
    if isempty(theta2)
        theta2 = pv.theta2
    end
    imp = if bnd == open
        pv_imp(pv.obs[2:end-1], [t.obs[2:end-1] for t in pt]..., cv=cv,L=L, theta1=theta1, theta2=theta2, bnd =open)
    elseif bnd ==periodic
        pv_imp(pv.obs, [t.obs for t in pt]..., cv=cv,L=L, theta1=theta1, theta2=theta2, bnd =periodic)
    end
    return juobs.Corr(imp,pv.kappa, pv.mu,pv.gamma,pv.y0,pv.theta1,pv.theta2,)
end

@doc raw"""
     pv0_imp(pv0, pt...; theta1,theta2, cv,L::Int64=1, bnd::Boundary=open)
     pv0_imp(pv0::juobs.Corr,pt::juobs.Corr ...;cv,L=1,theta1 = Float64[], theta2 = Float64[], bnd::Boundary=open)::juobs.Corr

Improve the correlator G_PV0 and G_PV0 with the tensor correlator according to the improvement equation

```math
    V_0^I(t,\vec p) = V_0(t,\vec p) - c_V \sum_{k=1}^3 i sin(ap^k) T_{0k}(t,\vec p)
```


# Arguments
  - pvo: correlator G_{PV0}
  - pt: correlators G_{PT}. It contains the correlator G_{PT01}, G_{PT02}, G_{PT03}

# Keyword Arguments
  - `cv`: mandatory keyword argument. It is the improvement coefficient
  - `L::Int64=1`: Lattice size. Needed to convert the theta angle into momenta.
  - `theta1=Float64[]`: theta angle of the first quark. If empty, is read from `pv`
  - `theta2=Float64[]`: theta angles of the second quark. If empty, is read from `pv`
  - `bnd` : Boundary condition of the lattice. Used to compute the derivatives

See also [`sym_der`](@ref), [`Boundary`](@ref)
"""
function pv0_imp(pv0, pt...; theta1,theta2, cv,L::Int64=1, bnd::Boundary=open)
    p = (theta1.-theta2)./L

    if all(p.==0)
        return pv0
    end

    if length(pt) !=3
        error("[pv0_imp] unexpected number of tensor current correlator")
    end

    dsin = sin.(p)
    aux = sin(p[1])*pt[1]+sin(p[2])*pt[2]+sin(p[3])*pt[3]
    return pv0.-cv.*aux
end

function pv0_imp(pv0::juobs.Corr,pt::juobs.Corr ...;cv,L=1,theta1 = Float64[], theta2 = Float64[], bnd::Boundary=open)
    check_corr(pv0,pt...,no_gamma)
    if isempty(theta1)
        theta1 = pv0.theta1
    end
    if isempty(theta2)
        theta2 = pv0.theta2
    end
    imp = if bnd ==open
        pv0_imp(pv0.obs[2:end-1], [x.obs[2:end-1] for x in pt]..., cv=cv, L=L,theta1=theta1, theta2=theta2, bnd=bnd)
    elseif bnd ==periodic
        pv0_imp(pv0.obs, [x.obs for x in pt]..., cv=cv, L=L,theta1=theta1, theta2=theta2, bnd=bnd)
    end
    return juobs.Corr(imp, pv0.kappa,pv0.mu, pv0.gamma,pv0.y0,pv0.theta1,pv0.theta2)
end
