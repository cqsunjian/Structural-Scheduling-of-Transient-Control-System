function G=GROPNORM(W,ND,N)
%gradient of the group sparsity regulation
G=zeros(size(W));
if(N==1)
    for j=1:size(W,2)/ND
        col_s=(j-1)*ND+1;
        col_end=j*ND;
        G(:,col_s:col_end)=1/norm(norm(W(:,col_s:col_end)));
    end 
else
for i=1:N
    for j=1:N
        row_s=(i-1)*ND*N+1;
        row_end=i*ND*N;
        col_s=(i-1)*ND*N+(j-1)*ND+1;
        col_end=(i-1)*ND*N+j*ND;
        G(row_s:row_end,col_s:col_end)=1/norm(norm(W(row_s:row_end,col_s:col_end)));
    end
end
end
end