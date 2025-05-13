function [res]=Int_txinf(b,A,t1,f)

res=-b*A^(-2)*(A*t1-eye(length(A)))*f;



end