function [res]=Int_tx(b,A,t1,t2,f)

% res=b*A^(-2)*expm(A*(t2-t1))*(t2*A-eye(length(A)))*f-b*A^(-2)*(A*t1-eye(length(A)))*f;




res=b*A^(-2)*expm(A*(t2-t1))*(t2*A-eye(length(A)))*f-b*A^(-2)*(A*t1-eye(length(A)))*f;


end