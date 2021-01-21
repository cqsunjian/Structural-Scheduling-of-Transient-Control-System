%simulation of the reinforcement learning control method A2C implemented in the power system transient
%control, with the contingency happpend at 0.1s and cleared at 0.4s
clc;

tic;
mpc = runpf('case118');
n=size(mpc.bus,1);
m=size(mpc.gen,1);
b=size(mpc.branch,1);
param=cell(n,1);
Pg=zeros(m,1);
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
% L=[-0.2,-0.9;0.01, 1];
abs(eig(L))
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
u0=zeros(2*m,1);
x=x0;
u=u0;
q=zeros(m,1);
p=zeros(m);
next_x=x;
next_u=u;

%action-critic network basis dimension number
ND=20;
% weight of action nerual network
W_ai=cell(m,1);
Ma=cell(m,1);

% weight of critic nerual network
W_ci=cell(m,1);
Wc=cell(m,1);

load('data\initialparam_A2C.mat');

% the weights of x_i and u_i
a_x=0.97;
a_u=0.03;

% len = 1.5/0.01 = 150
len=floor(last_time/samping_time);

x_record=zeros(size(x,1),len);
u_record=zeros(size(u,1),len);
q_record=zeros(m,len);

a=0.8;
Na=8;

%Umax = 5
Umax=5;

%learning rates
alpha=0.11;
beta=0.00001;
th=0.5;

%reference bus: bus 30
rb=find(mpc.bus(:,2)==3);
rb=find(mpc.gen(:,1)==rb);
PY=Y;

%start simulation
for k=2:len

    %start contingency,start at 1s
    if(k==floor(fault_time/samping_time))
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
    
    %start controller
    if(k>control_time/samping_time)
        xp=zeros(size(x,1)+m,1);
        xp(1:4:end) = x(1:3:end);
        xp(2:4:end) = x(2:3:end);
        xp(3:4:end) = x(3:3:end);
        
        %Q_i is a reinforcement signal
        xp(4:4:end) = q_record(:,k-1);
        xp(2:4:end)=xp(2:4:end)-Va;
        phi=PHI_A2C(xp,c,ND);
        for j=1:m
            if j==rb
                continue;
            end
            % calculate the output of actor network
            u(2*j-1:2*j)=Umax*tanh(Ma{j}*W_ai{j}*phi/Umax);
            
            %calculate the output of critic network
            q(j)=Wc{j}*W_ci{j}*phi;
        end
    end
    
    %save the k-th state, action and Q
    x_record(:,k)=x;
    u_record(:,k)=u;
    q_record(:,k)=q;
    for j=1:m
        if(j==rb)
            continue;
        end
        p(j)=norm(x_record((j-1)*3+1:j*3,k-1))*a_x+norm(u_record(j,k-1))*a_u>=th;
    end
    
    % use the model of power system to calculate the next x: next_x
    k1=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u,Vref);
    next_x=x+samping_time*k1;

    %updating the weights : the weights of actor network and the weights of critic network 
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
            
            %updating the weights of actor network
            W_ai{j} = W_ai{j}-alpha*Ma{j}'*diag(1-tanh(Ma{j}*W_ai{j}*phi/Umax).^2)*[1,0,0;0,0,1]*(next_y-L*y)*phi';
            
            %updating the weights of critic network 
            W_ci{j} =W_ci{j} - beta*(Wc{j}'*q(j)+p(j)*a^Na)*phi';
        end
    end
    
    x=next_x;
    
end
t=0:samping_time:samping_time*len-samping_time;

subplot(3,2,1:2);plot(t,x_record(1:3:end,:)');ylabel('angular frequency(rad/s)');xlabel('time(s)');%axis([0 5,-1.5,1.5]);
subplot(3,2,3:4);plot(t,x_record(2:3:end,:)');ylabel('phase');xlabel('time(s)');
subplot(3,2,5:6);plot(t,x_record(3:3:end,:)');ylabel('voltage Eq');xlabel('time(s)');
toc;