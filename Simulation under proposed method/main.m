%simulation of the proposed method
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
last_time=1;
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



%action network basis dimension number
ND=20;
range=diag([1,1,1])*[-10,10;-2 2;-1,1];
load('data\\initialparam.mat','W','Ma','Wc','c');
G=cell(m,1);
for i=1:m
   a_x=0.97;
   a_u=0.03;
   G{i}=zeros(size(W{i}));
end
len=floor(last_time/samping_time);
x_record=zeros(size(x,1),len);
u_record=zeros(size(u,1),len);
q_record=zeros(m,len);
e_record=zeros(m,len);
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

    %start contingency
    if(k==floor(fault_time/samping_time))
        %choose branch and change admittance
        mpc.branch(11,11)=0;
        mpc.branch(12,11)=0;
        Y=ReducedY(mpc);
    end
    %end contingency,
    if(k==floor(fault_clear_time/samping_time))
        %restore Y 
        mpc.branch(11,11)=1;
        mpc.branch(12,11)=1;
        Y=PY;
    end
    %start controller
    if(k>control_time/samping_time)
        xp=x;
        xp(2:3:end)=xp(2:3:end)-Va;
        phi=PHI(xp,c,ND);

    for j=1:m
        if j==rb
            continue;
        end
         u(2*j-1:2*j)=Umax*tanh(Ma{j}*W{j}*phi/Umax);%calculate control output by action neural network
         q(j)=Wc{j}*W{j}*phi;
         G{j}=GROPNORM(W{j},ND,1);
    end
%     disp(u(1:10)')
    end
    %record datas
    x_record(:,k)=x;
    u_record(:,k)=u;
    for j=1:m
        if(j==rb)
            continue;
        end
        p(j)=norm(x_record((j-1)*3+1:j*3,k-1))*a_x+norm(u_record(j,k-1))*a_u>=th;
    end
    x_record(2:3:end,k)=x_record(2:3:end,k)-Va;
    %calculate next_x
    k1=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u,Vref);
    next_x=x+samping_time*k1;
    %update paramters / learning
    if(k>floor(control_time/samping_time))
        for j=1:m
            if(j==rb)
                continue;
            end
            y=x((j-1)*3+1:(j-1)*3+3);
            y(2)=y(2)-Va(j);
            y(3)=y(3)-Eqs(j);
            next_y=next_x((j-1)*3+1:(j-1)*3+3);
            next_y(2)= next_y(2)-Va(j);
            next_y(3)=next_y(3)-Eqs(j);
            W{j}=W{j}-alpha*Ma{j}'*diag(1-tanh(Ma{j}*W{j}*phi/Umax).^2)*[1,0,0;0,0,1]*(next_y-L*y)*phi'-beta*(Wc{j}'*q(j)+p(j)*a^Na)*phi'-gamma*W{j}.*G{j};
            w=0;
            for o=1:size(W{j},2)/ND
                col_s=(o-1)*ND+1;
                col_end=o*ND;
                w=w+norm(norm(W{j}(:,col_s:col_end)))^2;
            end 
            e_record(j,k)=alpha*norm(next_y-L*y)^2/2+beta*norm(Wc{j}'*q(j)+p(j)*a^Na)^2/2+gamma*w/2;
        end
        
    end
  
    x=next_x;%update state x
    
end
t=0:samping_time:samping_time*len-samping_time;

plotW;
figure;
plot(t,e_record);xlabel('time(s)');ylabel('objective function J_i');axis([0.4,1,0,0.2]);grid on;
figure;
subplot(2,1,1);plot(t,x_record(1:3:end,:)');ylabel('angular frequency \omega(rad/s)');xlabel('time(s)');axis([0.01 1,-3,3]);grid on;
subplot(2,1,2);plot(t,x_record(3:3:end,:)');ylabel('voltage amplitude Eq(p.u.)');xlabel('time(s)');axis([0.01 1,0.8,1.1]);grid on;
