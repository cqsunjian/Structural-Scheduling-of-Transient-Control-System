function L = CalculatingL(Y,Va,Vm,M)
%calculation for matrix A
m=size(Y,1);
L=zeros(2*m,2*m);
PHI=zeros(m,m);
for i=1:m
    for j=1:m
        if(i==j)
            continue;
        end
        L(2*(i-1)+1,2*(j-1)+2)=-abs(Y(i,j))*Vm(i)*Vm(j)*sin(Va(i)-Va(j)-angle(-Y(i,j)));
         L(2*(i-1)+1,2*(i-1)+2)= L(2*(i-1)+1,2*(i-1)+2)- L(2*(i-1)+1,2*(j-1)+2);
    end
end


end

