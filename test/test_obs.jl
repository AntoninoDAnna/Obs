using Revise, ObsIO, Obs, ADerrors

Base.zero(::Type{uwreal}) = uwreal(0.0)
Base.abs(u::uwreal) = sqrt(u^2)

## Open BC
pa0  = ObsIO._read_corr("data/H101_x0_40_G5_G0G5_kappa_0.13675962_0.13675962_mu_0.0_0.0_theta1_0.0_0.0_0.0_theta2_0.0_0.0_0.0.bdio");
pp   = ObsIO._read_corr("data/H101_x0_40_G5_G5_kappa_0.13675962_0.13675962_mu_0.0_0.0_theta1_0.0_0.0_0.0_theta2_0.0_0.0_0.0.bdio");
a0a0 = ObsIO._read_corr("data/H101_x0_40_G0G5_G0G5_kappa_0.13675962_0.13675962_mu_0.0_0.0_theta1_0.0_0.0_0.0_theta2_0.0_0.0_0.0.bdio");

Obs.meff(pp)
Obs.mpcac(pa0,pp)
