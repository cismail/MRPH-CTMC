function [tau,Age,Cost] = PT_bs(Q,mu,b,max)

if nargin<4
    max=10;
end

eps=0.001;
tau=0.01;
max_tau=max;
min_tau=0.01;
cond=true;
lambda=1/tau*ones(length(Q),1);


[Age,Cost] =CTMC_anl_serv(Q,lambda,mu);

if Cost<b
    cond=false;   
end



while cond
tau=min_tau/2+max_tau/2;
lambda=1/tau*ones(length(Q),1);

[Age,Cost] =CTMC_anl_serv(Q,lambda,mu);

if Cost<b
    max_tau=tau;
else
    min_tau=tau;
end

cond=abs(Cost-b)>eps;
if tau==max
    warning('max is reached')
    Cost=-1;
    cond=false;
end


end

end

