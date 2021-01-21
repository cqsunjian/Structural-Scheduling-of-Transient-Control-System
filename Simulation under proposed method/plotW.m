%drawing the control structure
N=54;
n=54;
WX=zeros(54,54);
ac=12;%maximal allowable communication links for each
pc=20;%maximal allowable sensors
ac=25;%maximal allowable actuators
for i=1:N
    for j=1:N
        WX(i,j)=(sum(sum(W{i}(1:3,(j-1)*ND+1:j*ND).^2)));
        if(WX(i,j)<0.001)
            WX(i,j)=0;
        end
    end
end
WX(:,30)=0;
WX(30,:)=0;
[inx,jnx]=findmaxNindex(WX,ac*n);
CV=zeros(N,N);
CV(inx,jnx)=1;
WX=WX.*CV;

inx=findmaxNrow(WX,ac);
CV=zeros(N,N);
CV(inx,:)=1;
WX=WX.*CV;
jnx=findmaxNcol(WX,pc);
CV=zeros(N,N);
CV(:,jnx)=1;
WX=WX.*CV;

% figure;
spy(WX);

xlabel('remote node index');
ylabel('local node index');
% h=zlabel('F norm of blocks $\hat{X}_{ji}$');
% set(h,'Interpreter','latex','fontsize',12);