using Obs

cP  = 0.1;
cA0 = -0.2;
E = 0.5;
pp= [cP*exp(-E*t) for t in 1:40];
pa0 = [cA0*exp(-E*t) for t in 1:40];

mpcac = cA0/(2cP) * sinh(E);
fps   =  abs(cA0)/sqrt(cP)*sqrt(2/E);

ca = 0.02;
pa0_ca = pa0 + ca* Obs.sym_der(pp,Obs.open);

#improvements
let
    Obs.pa0_imp(pa0_ca,pp,ca=ca);
    Obs.pv_imp(pa0_ca,pp,cv=ca, theta1 = zeros(3), theta2 = zeros(3));
    Obs.pv_imp(pa0_ca,pp,pp,pp,pp, cv=ca,theta1 = [2/3,0.0,1/3],theta2 =[0,-1,0])
    Obs.pv0_imp(pa0_ca,pp,cv=ca,theta1 = zeros(3), theta2=zero(3),bnd=Obs.open)
    Obs.pv0_imp(pa0_ca,pp,pp,pp,cv=ca,theta1 = fill(1/3,3), theta2=zero(3),bnd=Obs.open)
    Obs.a0a0_imp(pa0_ca,0.5*pp,ca=ca)
    Obs.v_imp(pa0_ca,pp,cv=ca,theta1=zeros(3),theta2=zeros(3))
    Obs.v_imp(pa0_ca,pp,pp,pp,pp,pp,pp,pp,cv=ca,theta1=fill(1/7,3),theta2=zeros(3))
end
let
    Obs.mpcac(pa0,pp)
    Obs.meff(pa0)
    Obs.mpcac(pa0_ca,pp,ca)
    Obs.ps_dec(pa0,pp,E,0)
    Obs.dec(pa0,E,0)
end

true;
