function Q = gram_schmidt(A)
% ****** Gram-Schmidt method*****************
% by Jamie
% 07/23/2014
% *******************************************
[m,n]    =  size(A)                         ;
Q        =  zeros(m,n)                      ;
R        =  zeros(n,n)                      ;
for col  =  1:size(A,2)                     ;
    v    =  A(:,col)                        ;
    
    for row    =  1:col-1
        R(row,col) =  Q(:,row)'*A(:,col)    ;
        v          =  v-R(row,col)*Q(:,row) ;
    end
    
    R(col,col) =  norm(v)                   ;
    Q(:,col)   =  v/R(col,col)              ;
end

end
% *******************************************