function [res]=Int_tx22inf(b,A,t,f)



res=b*(-A^(-3)*(A^2*t^2-2*A*t+2*eye(length(A)) ))*f/2;


end