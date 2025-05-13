function [res]=Int_x22(b,A,d,f)

res= b*(expm(A*d)*(A^2*d^2-2*A*d+2*eye(length(A)))*A^(-3)-2*A^(-3))*f/2;
end