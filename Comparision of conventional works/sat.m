function y=sat(x,Um)
%Satuaration function
y=x;
if(x>Um)
    y=Um;
elseif(x<-Um)
    y=-Um;
end
end