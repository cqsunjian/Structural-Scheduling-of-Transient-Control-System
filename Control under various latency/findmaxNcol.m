function [jn]=findmaxNcol(W,nx)
%find nx columns with the largest sums of columns 
[m,n]=size(W);
    wv=sum(W,1);
    [ws,v]=sort(wv,'descend');
     for i=1:nx
        jn(i)=v(i);
     end
end