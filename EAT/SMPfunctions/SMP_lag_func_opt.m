function [obj,Table,meanAoII,MeanSamplingRate] = SMP_lag_func_opt(Q,mu,lagr,opt)

eps_value_dist=opt.eps_value_dist;
eps_pol_dist=opt.eps_pol_dist;

N=length(Q);
initial_tau=opt.initial_pol;

% Initial Policy Evaluation
[r,P,h,pik,meanAoII,MeanSamplingRate] = reward_lambda_tau(Q,mu,initial_tau,lagr);
Phat=eye(N)-P;
% for k=1:N
%     Phat(k,k)=1-P(k,k);
% end
A=[h',Phat(:,1:N-1)];
res=A\r';
g=res(1);
v=[res(2:N);0];
old_tau=initial_tau;
cond=true;

Age=[meanAoII];
Cost=MeanSamplingRate;
Taus=old_tau;
gT=g;
vT=v';
AgeLambdaCost=meanAoII+lagr*MeanSamplingRate;
%tic
while cond
    % Policy Improvement
    for i=1:N
        taus=old_tau;
        [new_a,mv] = policy_solver(Q,mu,lagr,taus,v,i,opt);
        new_tau(i)=new_a;
    end
    % Policy Evaluation
    [r,P,h,pik,meanAoII,MeanSamplingRate] = reward_lambda_tau(Q,mu,new_tau,lagr);
    Phat=eye(N)-P;
    A=[h',Phat(:,1:N-1)];
    res=A\r';
    g=res(1);
    v=[res(2:N);0];
    
    pol_dist=norm(old_tau-new_tau);
    cond1 = (pol_dist > eps_pol_dist);
    cond2= (norm(gT(end)-g))>(eps_value_dist);
    cond3=size(Taus,1)<30;
    if not(cond3)
        warning('Maximum number of iteration is exceeded')
    end
    cond=and(and(cond1,cond2),cond3);
    old_tau=new_tau;

    Age=[Age;meanAoII];
    Cost=[Cost;MeanSamplingRate];
    Taus=[Taus;old_tau];
    gT=[gT;g];
    vT=[vT;v'];
    AgeLambdaCost=[AgeLambdaCost;meanAoII+lagr*MeanSamplingRate];

    Table=table(AgeLambdaCost,Age,Cost,Taus,gT,vT);
 %   disp(Table)

end
%toc
obj=meanAoII+lagr*MeanSamplingRate;
%disp(Table)

end

