function En = CalculateEnergy(Y,Va,Vac,Vm)
%energy calculation on branches
Yb=abs(imag(Y));
m=size(Y,1);
En=zeros(m,1);
for i=1:m
    if i==30
        continue;
    end
    for j=1:m
       En(i)=En(i)+Vm(i)*Vm(j)*Yb(i,j)*(cos(Va(i)-Va(j))-cos(Vac(i)-Vac(j))); 
    end
end
end

