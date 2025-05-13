function [res]=Int_1(b,A,d,f)

res=b*inv(A)*(expm(A*d)-eye(length(A)))*f;

end