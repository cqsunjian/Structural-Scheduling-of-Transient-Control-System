function y=bas(x,c,w,seplength)
%radial basis function 
l=size(c,1)/seplength;%subsystem number
n=size(c,2);%dimension of radias basis
y=zeros(l*n,1);%output
for i=1:l% for each subsystem
    for j=1:n %for eacg dimension of basis
        y(j+(i-1)*n)=exp(-norm(x((i-1)*seplength+1:i*seplength,1)-c((i-1)*seplength+1:i*seplength,j))/w^2);
    end
end
end