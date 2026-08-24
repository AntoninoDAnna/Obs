using Revise, ObsIO, Obs, ADerrors

Base.zero(::Type{uwreal}) = uwreal(0.0)

## Open BC
pa0  = ObsIO._read_corr("data/H101_x0_40_G5_G0G5_kappa_0.13675962_0.13675962_mu_0.0_0.0_theta1_0.0_0.0_0.0_theta2_0.0_0.0_0.0.bdio");
pp   = ObsIO._read_corr("data/H101_x0_40_G5_G5_kappa_0.13675962_0.13675962_mu_0.0_0.0_theta1_0.0_0.0_0.0_theta2_0.0_0.0_0.0.bdio");
a0a0 = ObsIO._read_corr("data/H101_x0_40_G0G5_G0G5_kappa_0.13675962_0.13675962_mu_0.0_0.0_theta1_0.0_0.0_0.0_theta2_0.0_0.0_0.0.bdio");

# improve pa0
let
    c = Obs.pa0_imp(pa0,pp,ca=1)
    der = [uwreal(0.0);Obs.sym_der(pp.obs[2:end-1],Open);uwreal(0.0)]
    ci =pa0.obs .- der
    aux  = c.obs .- ci
    uwerr.(aux)
    if any(x->abs(x.mean)>eps(), aux)
        error("error in impove pa0 with open  open boundary conditions")
    end
end

# improve a0a0
let
    c = Obs.a0a0_imp(a0a0,pa0,ca=1)
    der = [uwreal(0.0);Obs.sym_der(pa0.obs[2:end-1],Open);uwreal(0.0)]
    ci =a0a0.obs .- 2der
    aux  = c.obs .- ci
    uwerr.(aux)
    if any(x->abs(x.mean)>eps(), aux)
        error("error in impove a0a0 with open  open boundary conditions")
    end
end


# improve vivi
let
    c = Obs.v_imp(a0a0,pa0,cv=1)
    der = [uwreal(0.0);Obs.sym_der(pa0.obs[2:end-1],Open);uwreal(0.0)]
    ci =a0a0.obs .- 2der
    aux  = c.obs .- ci
    uwerr.(aux)
    if any(x->abs(x.mean)>eps(), aux)
        error("error in impove vivi with open  open boundary conditions")
    end
end

# improve pvi
let
    c = Obs.pv_imp(pa0,pp,cv=1,real=true)
    der = [uwreal(0.0);Obs.sym_der(pp.obs[2:end-1],Open);uwreal(0.0)]
    ci =pa0.obs .- 2der
    aux  = c.obs .- ci
    uwerr.(aux)
    if !any(x->abs(x.mean)>eps(), aux)
        error("error in impove pvi with open  open boundary conditions")
    end
end

# improve pv0
let
    c = Obs.pv0_imp(pa0,pp,pp,pp,cv=1,real=true)
    der = [uwreal(0.0);Obs.sym_der(pp.obs[2:end-1],Open);uwreal(0.0)]
    ci =pa0.obs .- 2der
    aux  = c.obs .- ci
    uwerr.(aux)
    if !any(x->abs(x.mean)>eps(), aux)
        error("error in impove pv0 with open  open boundary conditions")
    end
end
