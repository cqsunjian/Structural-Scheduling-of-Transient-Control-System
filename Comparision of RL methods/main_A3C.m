%simulation of the reinforcement learning control method A3C implemented in the power system transient
%control, with the contingency happpend at 0.1s and cleared at 0.4s
clear;
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
p=zeros(m,1);


load('data\initialparam_A3C.mat');
len=floor(last_time/samping_time);
x_record=zeros(size(x,1),len);
u_record=zeros(size(u,1),len);

a=0.8;
Na=8;
Umax=5;
ND=20;
alpha=0.09;
beta=0.00001;
gamma=0.001; 
a_x=0.97;
a_u=0.03;
l=0.1;
th=0.5;%cost threshold
op =0;

rb=find(mpc.bus(:,2)==3);
rb=find(mpc.gen(:,1)==rb);

PY=Y;

%Initialize 4 threads the first threads is main threads
number = 4;

% Initialize the size of parameters in each thread
u_par=cell(number,1);
p_par = cell(number,1);
q_par = cell(number,1);
phi = cell(number,1);
W_ai_par = cell(number,1);
W_ci_par = cell(number,1);
Ma_par = cell(number,1);
Wc_par = cell(number,1);
xp = cell(number,1);
u_record_par = cell(number,1);
q_record_par = cell(number,1);
k_par=zeros(size(x,1),number);
next_x=zeros(size(x,1),number);

% Initialize the parameters in each thread
for i =1:number
    W_ai_par{i} =cell(m,1);
    W_ci_par{i} =cell(m,1);
    Ma_par{i}    = cell(m,1);
    Wc_par{i}    = cell(m,1);
    for j=1:m
        W_ai_par{i}{j} = (rand(3,ND*m)-0.5)/(3*ND*m*ND*m); 
        W_ci_par{i}{j} = (rand(3,ND*m)-0.5)/(3*ND*m*ND*m);
        Wc_par{i}{j}    = (rand(1,3)-0.5);
        Ma_par{i}{j}    = (rand(2,3)-0.5);
    end
    u_par{i} = zeros(size(u,1),1);
    phi{i} = zeros(m*ND,1);
    p_par{i} =zeros(m,1);
    q_par{i} = zeros(m,1);
    xp{i} = zeros(size(x,1),1); 
    u_record_par{i} = zeros(size(u,1),len);
    q_record_par{i} = zeros(m,len);
end

%start simulation
for k=2:len

    %start contingency
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
    %start controller
    if(k>control_time/samping_time)
        % Start 4 threads
        parfor num =1:number
            xp{num}=zeros(size(x,1)+m,1);
            xp{num}(1:4:end) = x(1:3:end);
            xp{num}(2:4:end) = x(2:3:end);
            xp{num}(3:4:end) = x(3:3:end);
            xp{num}(4:4:end) = q_record_par{num}(:,k-1); %Q_i is a reinforcement signal
            xp{num}(2:4:end)=xp{num}(2:4:end)-Va;
            phi{num}=PHI_A3C(xp{num},c,ND);
            for j=1:m
                if j==rb
                    continue;
                end
                %Calculate the output of actor network in each thread
                u_par{num}(2*j-1:2*j)=Umax*tanh(Ma_par{num}{j}*W_ai_par{num}{j}*phi{num}/Umax);
                
                %Calculate the output of critic network in each thread
                q_par{num}(j)=Wc_par{num}{j}*W_ci_par{num}{j}*phi{num};
            end 
        end
    end
   
    %record datas
    x_record(:,k)=x;
    
    for num=1:number
        u_record_par{num}(:,k)=u_par{num};
        q_record_par{num}(:,k) = q_par{num};
        for j=1:m
            if(j==rb)
                continue;
            end
            p_par{num}(j)=norm(x_record((j-1)*3+1:j*3,k-1))*a_x+norm(u_record_par{num}(j,k-1))*a_u>=th;
        end
        k_par(:,num)=nonlinear_dynamic(D,M,mpc,GenBus,Y,x,u_par{num},Vref);
        next_x(:,num)=x+samping_time*k_par(:,num);
    end
    
    sum =zeros(number,1);
    for num=1:number
        for j=1:m
            if(j==rb)
                continue;
            end
            sum(num) =sum(num) + norm(next_x((j-1)*3+1:j*3,num));
        end
    end
    
    flag =1;
    Min =sum(1);
    for num =2:number
        if Min>=sum(num)
            Min = sum(num);
            flag = num;
        end
    end
    u_record(:,k) = u_par{flag};
    
    %update paramters / learning
    if(k>floor(control_time/samping_time))
        for num =1:number
            for j=1:m
                if(j==rb)
                    continue;
                end
                y=x((j-1)*3+1:(j-1)*3+3);
                y(2)=y(2)-Va(j);
                y(3)=y(3)-Eqs(j);
                next_y=next_x((j-1)*3+1:(j-1)*3+3,flag);
                next_y(2)= next_y(2)-Va(j);
                next_y(3)=next_y(3)-Eqs(j);
                % Update the weights of actor network in each thread
                W_ai_par{num}{j}=W_ai_par{num}{j} - alpha*Ma_par{num}{j}'*diag(1-tanh(Ma_par{num}{j}*W_ai_par{num}{j}*phi{num}/Umax).^2)*[1,0,0;0,0,1]*(next_y-L*y)*phi{num}';
                
                % Update the weights of critic network in each thread
                W_ci_par{num}{j} =W_ci_par{num}{j} - beta*(Wc_par{num}{j}'*q_par{num}(j)+p_par{num}(j)*a^Na)*phi{num}';
           end
        end
    end
    %choose the best weight
    W_ai_par{1} = W_ai_par{flag};
    W_ci_par{1} = W_ci_par{flag};    
    x=next_x(:,flag);
    
end
t=0:samping_time:samping_time*len-samping_time;
subplot(3,2,1:2);plot(t,x_record(1:3:end,:)');ylabel('angular frequency(rad/s)');xlabel('time(s)');%axis([0 5,-1.5,1.5]);
subplot(3,2,3:4);plot(t,x_record(2:3:end,:)');ylabel('phase');xlabel('time(s)');
subplot(3,2,5:6);plot(t,x_record(3:3:end,:)');ylabel('voltage Eq');xlabel('time(s)');
toc;