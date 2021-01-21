function phi=PHI_A2C(x,c,ND)
%radiaul basis function for A2C
m=size(x,1)/4;
phi=zeros(ND*m,1);
sigma=20;
for i=1:m
    xi=x((i-1)*4+1:(i-1)*4+4); 
    for j=1:ND
        cj=c(:,j);
        phi((i-1)*ND+j)=exp(-norm(xi-cj)^2/sigma/2);
    end
end  
phi=phi/sqrt(m/ND);;
end