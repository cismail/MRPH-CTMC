function [res]=Int_1inf(b,A,f)

res=b*inv(A)*(-eye(length(A)))*f;

end