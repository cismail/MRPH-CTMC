function [tau,Age,Cost] = ST_bs(Q,mu,b,max)

if nargin<4
    max=10;
end

eps=0.001;
tau=0;
max_tau=max;
min_tau=0;
cond=true;
%[reward,P,cycleVec,pik,Age,Cost] = reward_lambda_tau(Q,mu,tau,0);

[Age,Cost] =CTMC_anl_tau(Q,mu,tau);

if Cost<b
    cond=false;   
end



while cond
tau=min_tau/2+max_tau/2;

[Age,Cost] =CTMC_anl_tau(Q,mu,tau);

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

