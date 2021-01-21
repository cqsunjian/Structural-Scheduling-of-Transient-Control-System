function dx=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u,Vref)
%power system dynamic model, where the output is the time derivative of
%state x
rb=find(mpc.bus(:,2)==3);
rb=find(mpc.gen(:,1)==rb);
dx=zeros(size(x));
xd=0.138;
xq=0.0396;
xd2=0.0396;%xd'
Xd=xd-xd2;
Kf=10;%gain of AVR
Td=3;
% Eqs=mpc.bus(GenBus,8);%steady state values of Edi???



Pd=mpc.bus(GenBus,3)/mpc.baseMVA;
Pg=mpc.gen(:,2)/mpc.baseMVA;
% n=size(x,1)/3;
% Eqs=mpc.bus(GenBus,8);%steady state values of Edi???
Eqv=(x(3:3:end)).*exp(1j*x(2:3:end));%calculate internal voltage??
% Xd=ones(n,1);%vector of xdi-xdi'
% Y=imag(Y);
S=Eqv.*conj(Y*Eqv);
PE=real(S);
QE=imag(S);
Eq=x(3:3:end);



vq=Eq-xd*(u(2:2:end)-QE)./Eq;%q-voltage
vd=xq*(u(1:2:end)-PE)./Eq;%d-voltage
Vt=sqrt(vq.^2+vd.^2);%terminal voltage

Ef=Kf*(Vref-Vt);%excitment voltage

dx(2:3:end)=x(1:3:end);
dx(1:3:end)=dx(1:3:end)+1/M*(-D*x(1:3:end)+Pg-Pd-PE+u(1:2:end));%+0.0*(rand(size(x)/2)-0.5); frequency dynamics
dx(3:3:end)=dx(3:3:end)+1/Td*(-Eq-Xd*(QE-u(2:2:end))./Eq+Ef);% iternal voltage dynamics
dx((rb-1)*3+1:(rb-1)*3+2)=0;%reference bus set

end