function [jn]=findmaxNcol(W,nx)
%find N column with largest sum of columns, with the reture is the column index 
[m,n]=size(W);
    wv=sum(W,1);
    [ws,v]=sort(wv,'descend');
     for i=1:nx
        jn(i)=v(i);
     end
end