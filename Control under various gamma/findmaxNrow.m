function [in]=findmaxNrow(W,nx)
%find nx rows with largest sum of rows. the return variabe is the row index
[m,n]=size(W);
    wv=sum(W,2);
    [ws,v]=sort(wv,'descend');
     for i=1:nx
        in(i)=v(i);
     end
end