function [in,jn]=findmaxNindex(W,nx)
% find nx elements with the largest values
[m,n]=size(W);
    wv=reshape(W,1,m*n);
    [ws,v]=sort(wv,'descend')
     for i=1:nx
         jn(i)=fix((v(i)-1)/m)+1;
         in(i)=mod(v(i)-1,m)+1;
     end
end