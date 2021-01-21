function phi=PHI(x,c,ND)
%radius basis function
m=size(x,1)/3;
phi=zeros(ND*m,1);
sigma=20;
for i=1:m
    xi=x((i-1)*3+1:(i-1)*3+3); 
    for j=1:ND
        cj=c(:,j);
        phi((i-1)*ND+j)=exp(-norm(xi-cj)^2/sigma/2);
    end
end  
phi=phi/sqrt(m/ND);
% phi=x;
% phi=phi;
%phi=repmat(phi,N,1);
end