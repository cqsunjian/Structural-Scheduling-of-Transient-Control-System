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
control_time=2;
last_time=1.5;
fault_time=0.1;
fault_clear_time=0.2;
samping_time=0.01;
D=0.002;
M=0.050;
dw=2.5;
xd=0.138;
xq=0.0396;
xd2=0.0396;%xd'
Xd=xd-xd2;
dt=0.2;
Td=3;
Kf=10;%gain of AVR
de=100;
L=[1-samping_time*(D+dw)/M,-samping_time*dt/M 0;...
    samping_time,1,0;...
    0,0,1-samping_time*(1+de)/Td];
[Y,GenBus]=ReducedY(mpc);
as=mpc.bus(GenBus,9)/180*pi;


Eqs=mpc.bus(GenBus,8);%steady state values of Eq
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
beta=0.001;
gamma=0.0001;
l=0.1;
th=0.5;%cost threshold
rb=find(mpc.bus(:,2)==3);
rb=find(mpc.gen(:,1)==rb);

PY=Y;

%start simulation
for k=2:len

    %start contingency,start at 1s
    if(k==floor(fault_time/samping_time))
        %generator bus
        mpc.branch(11,11)=0;
        mpc.branch(12,11)=0;
        Y=ReducedY(mpc);
    end
    %fault clear
    if(k==floor(fault_clear_time/samping_time))
        mpc.branch(11,11)=1;
         mpc.branch(12,11)=1;
        Y=PY;
    end
    %record datas
    x_record(:,k)=x;
    %calculate next_x
    k1=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u,Vref);
    next_x=x+samping_time*k1;
    x=next_x;%update state x  
end
t=0:samping_time:samping_time*len-samping_time;
subplot(3,2,1:2);plot(t,x_record(1:3:end,:)');ylabel('angular frequency \omega(rad/s)');xlabel('time(s)');axis([0.01 1.2,-50,50]);grid on;
subplot(3,2,3:4);plot(t,x_record(2:3:end,:)');ylabel('phase \theta(rad)');xlabel('time(s)');axis([0.01 1.2,-2,2]);grid on;
subplot(3,2,5:6);plot(t,x_record(3:3:end,:)');ylabel('voltage amplitude Eq(p.u.)');xlabel('time(s)');axis([0.01 1.2,0,1.1]);grid on;

