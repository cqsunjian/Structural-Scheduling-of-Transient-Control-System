function [in,jn]=findmaxNindex(W,nx)
%find N largest elements. the return variables are the row indexs and
%column indexs
[m,n]=size(W);
    wv=reshape(W,1,m*n);
    [ws,v]=sort(wv,'descend')
     for i=1:nx
         jn(i)=fix((v(i)-1)/m)+1;
         in(i)=mod(v(i)-1,m)+1;
     end
end