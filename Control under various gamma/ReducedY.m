function [Y,GenBus]=ReducedY(mpc)
%function for calculating admittance matrix 
GenBus=mpc.gen(:,1);
LodBus=setdiff(mpc.bus(:,1),GenBus);
[Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
VLm=mpc.bus(LodBus,8);
SL=(mpc.bus(LodBus,3)+1j*mpc.bus(LodBus,4))/mpc.baseMVA;
YL=conj(SL./VLm.^2);
m=size(Ybus,1);
Ybus=Ybus+sparse(LodBus,LodBus,YL,m,m);
%start Kron Reduction
Ynn=Ybus(GenBus,GenBus);
Ynr=Ybus(GenBus,LodBus);
Yrn=Ybus(LodBus,GenBus);
Yrr=Ybus(LodBus,LodBus);
Y=Ynn-Ynr*Yrr^-1*Yrn;
end