%Simulation of a decentralized droop control strategy with an energy function based ESS allocation scheme
clear;
clc;

mpc = runpf('case118');
n=size(mpc.bus,1);
m=size(mpc.gen,1);
b=size(mpc.branch,1);
param=cell(n,1);
Pg=zeros(m,1);
Vg=zeros(m,1);
theta=zeros(m,1);
control_time=0.4;
last_time=1.5;
fault_time=0.1;
fault_clear_time=0.2;
samping_time=0.01;
D=0.002;
M=0.050;
dw=5;
xd=0.138;
xq=0.0396;
xd2=0.0396;%xd'
Xd=xd-xd2;
dt=0.2;
Td=3;
Kf=10;%gain of AVR
de=100;
L=[1-samping_time*(D+dw)/M,-samping_time*dt/M;samping_time,1];

[Y,GenBus]=ReducedY(mpc);
as=mpc.bus(GenBus,9)/180*pi;

%set x0 to an SEP
Eqs=mpc.bus(GenBus,8);%steady state values of Eq
Vm=Eqs;
Va=mpc.bus(GenBus,9)/180*pi;
Vg=mpc.bus(GenBus,8).*exp(1j*Va);%internal voltage
Se=Vg.*conj(Y*Vg);
PEs=real(Se);
QEs=imag(Se);
vqs=Eqs-xd*(-QEs)./Eqs;%q-voltage
vds=xq*(-PEs)./Eqs;%d-voltage
vts=sqrt(vqs.^2+vds.^2);
Vref=(Kf*vts+Eqs+Xd*QEs./Eqs)/Kf;%terminal voltage

x0=zeros(3*m,1);
x0(2:3:end)=mpc.bus(GenBus,9)/180*pi;
Va=mpc.bus(GenBus,9)/180*pi;
x0(3:3:end)=Eqs;
%x0(2:2:end)=theta;
u0=zeros(2*m,1);

x=x0;
u=u0;
q=zeros(m);
p=zeros(m);
next_x=x;
next_u=u;



%action network basis dimension number
ND=5;
range=diag([1 0.1]*100)*[-1,1;-1 1];
W=cell(m,1);
Ma=cell(m,1);
Wc=cell(m,1);
Q=cell(m,1);
R=cell(m,1);
G=cell(m,1);
for i=1:m
   W{i}=(rand(ND*m,ND*m)-0.5)/ND/m; 
   Ma{i}=(rand(1,ND*m)-0.5);
   Q{i}=eye(2);
   R{i}=eye(1,1)*0.01;
   Wc{i}=(rand(1,ND*m)-0.5);
   G{i}=zeros(size(W{i}));
end

len=floor(last_time/samping_time);
x_record=zeros(size(x,1),len);
u_record=zeros(size(u,1),len);
q_record=zeros(m,len);
a=0.8;Na=8;
RECORDWEVA=zeros(len,1);
RECORDALPHAEVA=zeros(len,1);
RECORDPEVA=zeros(len,1);
Umax=5;
alpha=0.1;
beta=0.0001;
gamma=0.001;
l=0.1;
th=0.5;%threshold
rb=find(mpc.bus(:,2)==3);
rb=find(mpc.gen(:,1)==rb);

PY=Y;
kc=13;
Uon=2:28;
Vac=0;

%start simulation
for k=2:len

    %start contingency,start at fault time
    if(k==floor(fault_time/samping_time))
        %generator bus
        mpc.branch(11,11)=0;
        mpc.branch(12,11)=0;
        Y=ReducedY(mpc);
    end
    %fault clear
    if(k==floor(fault_clear_time/samping_time))
        %restore Y and cut off 26-25
        mpc.branch(11,11)=1;
         mpc.branch(12,11)=1;
        Y=PY;
        Vac=x(2:3:end);
        En=CalculateEnergy(Y,Va,Vac,Vm);
        [ens,ind]=sort(En,'descend');
        Uon=ind(1:20);
%         x0=rand(m,1);
    end
    %start controller
    if(k>control_time/samping_time)
        xp=x;
        xp(2:3:end)=xp(2:3:end)-Va;
    for j=1:m
        if j==rb
            continue;
        end
        if(ismember(j,Uon))
            u(2*(j-1)+1)=sat(-kc*x((j-1)*3+1),Umax);
        end
    end
    end
    %record datas
    x_record(:,k)=x;
    u_record(:,k)=u;
    k1=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u,Vref);
    next_x=x+samping_time*k1; 
    x=next_x;%update state x
    
end
t=0:samping_time:samping_time*len-samping_time;

subplot(3,2,1:2);plot(t,x_record(1:3:end,:)');ylabel('angular frequency \omega(rad/s)');xlabel('time(s)');%axis([0.01 1,-2.5,2.5]);grid on;
subplot(3,2,3:4);plot(t,x_record(2:3:end,:)');ylabel('\theta-\theta_s(rad)');xlabel('time(s)');%axis([0.01 1,-0.5,0.5]);grid on;
subplot(3,2,5:6);plot(t,x_record(3:3:end,:)');ylabel('voltage amplitude Eq(p.u.)');xlabel('time(s)');%axis([0.01 1,0.8,1.1]);grid on;
