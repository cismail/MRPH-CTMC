function [res]=Int_x(b,A,d,f)

res=b*A^(-2)*expm(A*d)*(d*A-eye(length(A)))*f+b*A^(-2)*f;

end