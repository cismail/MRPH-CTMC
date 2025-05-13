function [obj,Table,meanAoII,MeanSamplingRate] = SMP_lag_func(Q,mu,tmax,dt,lagr)


N=length(Q);
initial_tau=0*ones(1,N);

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
    a_space=[0:dt:tmax];
    for i=1:N
        taus=old_tau;
        [new_a,mv] = policy_improve_bisection(Q,mu,taus,lagr,v,a_space,dt,i);
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
    cond1 = (pol_dist > (N-1)*dt);
    cond2= (norm(gT(end)-g))>(0.001);
    cond=and(cond1,cond2);
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

