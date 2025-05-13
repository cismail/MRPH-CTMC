function [res]=Int_tx22(b,A,t1,t2,f)


term_1=expm(A*(t2-t1))*(A^2*t2^2-2*A*t2+2*eye(length(A)))*A^(-3);

term_2=(A^2*t1^2-2*A*t1+2*eye(length(A)))*A^(-3);

res=b*(term_1-term_2)*f/2;

end