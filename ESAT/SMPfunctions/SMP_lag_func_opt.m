function [obj,Table,meanAoII,MeanSamplingRate] = SMP_lag_func_opt(Q,mu,lagr,opt)
tvec=ones(3,1);
eps_value_dist=opt.eps_value_dist;
eps_pol_dist=opt.eps_pol_dist;

N=length(Q);
initial_tau=opt.initial_pol;

[meanAoII,MeanSamplingRate,r,h,P]=M3RPhaseType(Q,mu,initial_tau,lagr);
% Initial Policy Evaluation
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

Age=[tvec*meanAoII];
Cost=tvec*MeanSamplingRate;
Taus=old_tau;
gT=tvec*g;
vT=tvec*v';
AgeLambdaCost=tvec*(meanAoII+lagr*MeanSamplingRate);
%tic
while cond
    % Policy Improvement
    for i=1:N
        [new_a] = policyiter_es(Q,mu,lagr,i,v,old_tau);
        new_tau(i,:)=new_a;
    end
    % Policy Evaluation
    [meanAoII,MeanSamplingRate,r,h,P]=M3RPhaseType(Q,mu,new_tau,lagr);
    Phat=eye(N)-P;
    A=[h',Phat(:,1:N-1)];
    res=A\r';
    g=res(1);
    v=[res(2:N);0];
    
    pol_dist=norm(old_tau-new_tau);
    cond1 = (pol_dist > eps_pol_dist);
    cond2= (norm(gT(end)-g))>(eps_value_dist);
    cond=and(cond1,cond2);
    old_tau=new_tau;

    Age=[Age;tvec*meanAoII];
    Cost=[Cost;tvec*MeanSamplingRate];
    Taus=[Taus;old_tau];
    gT=[gT;tvec*g];
    vT=[vT;tvec*v'];
    AgeLambdaCost=[AgeLambdaCost;tvec*(meanAoII+lagr*MeanSamplingRate)];

    Table=table(AgeLambdaCost,Age,Cost,Taus,gT,vT);
    disp(Table)

end
%toc
obj=meanAoII+lagr*MeanSamplingRate;
%disp(Table)

end

