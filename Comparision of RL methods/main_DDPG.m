%simulation of the reinforcement learning control method DDPG implemented in the power system transient
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
last_time=1.0;
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
q=zeros(m,1);
p=zeros(m,1);
next_x=x;
next_u=u;

%action network basis dimension number
ND=20;
% weight of action nerual network
W_ai=cell(m,1);
Ma=cell(m,1);
G_ai=cell(m,1);

% weight of critic nerual network
W_ci=cell(m,1);
Wc=cell(m,1);
G_ci=cell(m,1);
load('Data\initialparam_DDPG.mat');

a_x=0.97;
a_u=0.03;

len=floor(last_time/samping_time);

%  epoches
epoches =3;

%save x£¬u£¬q
x_record=zeros(size(x,1),len);
u_record=zeros(size(u,1),len);
q_record=zeros(m,len);
p_record=zeros(m,len);

%arfa=0.8
a=0.8;

Na=8;

%Umax = 5
Umax=5;

%learning rates
alpha=0.13;
beta=0.00018;
gamma=0.001; 

tao = 0.001;

th=0.5;%cost threshold

%reference  bus = 30
rb=find(mpc.bus(:,2)==3);
rb=find(mpc.gen(:,1)==rb);
PY=Y;

%% start simulation
for k=2:len

    %start contingency,start at 1s
    if(k==floor(fault_time/samping_time))
        %choose branch and change admittance 26-25 near 26, they all are
        %generator bus
        mpc.branch(11,11)=0;
        mpc.branch(12,11)=0;
        Y=ReducedY(mpc);
    end
    
    %end contingency, end at 1.1s,and move to another SEP  
    if(k==floor(fault_clear_time/samping_time))
        %restore Y and cut off 26-25
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
        xp(4:4:end) = q_record(:,k-1);
        xp(2:4:end)=xp(2:4:end)-Va;
        phi=PHI_A2C(xp,c,ND);
        for j=1:m
            if j==rb
                continue;
            end
            % calculte the output of target_actor neural network
            u(2*j-1:2*j)=Umax*tanh(Ma{j}*W_ai{j}*phi/Umax);
            
            %calculate the output of target_critic neural network
            q(j)=Wc{j}*W_ci{j}*phi;
        end
    end
    %record datas
    x_record(:,k)=x;
    u_record(:,k)=u;
    q_record(:,k)=q;
    for j=1:m
        if(j==rb)
            continue;
        end
        p(j)=norm(x_record((j-1)*3+1:j*3,k-1))*a_x+norm(u_record(j,k-1))*a_u>=th;
    end
    p_record(:,k) = p;
    %calculate next_x
    k1=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u,Vref);
    next_x=x+samping_time*k1;

    %update paramters / learning
    
    W_ai_local = W_ai;
    W_ci_local = W_ci;
    x_local = x;
    u_local = u;
    q_local = q;
    p_local = p;
    x_record_local= zeros(size(x_local,1),epoches);
    u_record_local = zeros(size(u_local,1),epoches);
    q_record_local = zeros(size(q_local,1),epoches);
    p_record_local = zeros(size(p_local,1),epoches);
    
    for epoche = 2:epoches
        xp_local=zeros(size(x_local+m,1)+m,1);
        xp_local(1:4:end) = x_local(1:3:end);
        xp_local(2:4:end) = x_local(2:3:end);
        xp_local(3:4:end) = x_local(3:3:end);
        xp_local(4:4:end) = q_record_local(:,epoche-1);
        xp_local(2:4:end)=xp_local(2:4:end)-Va;
        phi=PHI_A2C(xp_local,c,ND);
        for num=1:m
            if num==rb
                continue;
            end
            % calculate the output of local_actor neural network 
            u_local(2*num-1:2*num)=Umax*tanh(Ma{num}*W_ai_local{num}*phi/Umax);
            
            % calculate the output of local_critic neural network 
            q_local(num)=Wc{num}*W_ci_local{num}*phi;
        end        
        x_record_local(:,epoche) = x_local;
        u_record_local(:,epoche) = u_local;
        q_record_local(:,epoche) = q_local;
        for num=1:m
            if(num==rb)
                continue;
            end
            p_local(num)=norm(x_record_local((num-1)*3+1:num*3,epoche))*a_x+norm(u_record_local(num,epoche))*a_u>=th;
        end
         p_record_local(:,epoche) = p;
         k1_local=nonlinear_dynamic(D,M,mpc,GenBus,Y,x_local,u_local,Vref);
         next_x_local=x_local+samping_time*k1_local;
         for j=1:m
            if(j==rb)
                continue;
            end
            y=x_local((j-1)*3+1:(j-1)*3+3);
            y(2)=y(2)-Va(j);
            y(3)=y(3)-Eqs(j);
            next_y=next_x_local((j-1)*3+1:(j-1)*3+3);
            next_y(2)= next_y(2)-Va(j);
            next_y(3)=next_y(3)-Eqs(j);   
            
            %update the weight of local_acator neural netwok 
            W_ai_local{j} = W_ai_local{j}-alpha*Ma{j}'*diag(1-tanh(Ma{j}*W_ai_local{j}*phi/Umax).^2)*[1,0,0;0,0,1]*(next_y-L*y)*phi';
            
            %update the weight of local_critic neural network
            W_ci_local{j} =W_ci_local{j} - beta*(a^(-1)*Wc{j}'*q_local(j)+p_local(j)*a^(Na+1)-q_record(j,epoche))*phi';
         end         
         
    end
    
    if(k>floor(control_time/samping_time))
        for j=1:m
            if(j==rb)
                continue;
            end
      
            %update the weight of target_acator neural netwok 
            W_ai{j} = tao*W_ai{j} + (1-tao)*W_ai_local{j};
            %update the weight of target_acator neural netwok 
            W_ci{j} = tao*W_ci{j} + (1-tao)*W_ci_local{j};
        end
    end
    
    x=next_x;%update state x
    
end
t=0:samping_time:samping_time*len-samping_time;

subplot(3,2,1:2);plot(t,x_record(1:3:end,:)');ylabel('angular frequency(rad/s)');xlabel('time(s)');%axis([0 5,-1.5,1.5]);
subplot(3,2,3:4);plot(t,x_record(2:3:end,:)');ylabel('phase');xlabel('time(s)');
subplot(3,2,5:6);plot(t,x_record(3:3:end,:)');ylabel('voltage Eq');xlabel('time(s)');
toc;