function [new_tau,mv] = policy_improve_bisection(Q,mu,taus,lagr,v,a,d,ind)

min_a=min(a);
max_a=max(a);
test_a=min_a/2+max_a/2;

taus(ind)=test_a;
[r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
test_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);

cond=true;

while cond
taus(ind)=test_a-d;
[r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
left_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);

if left_val<=test_val
    max_a=test_a;
    test_a=min_a/2+max_a/2;

    if abs(test_a-min(a))<d
        test_a=min(a);
        cond=false;
    end

    taus(ind)=test_a;
    [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
    test_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);

else 
taus(ind)=test_a+d;
[r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
right_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);
    if right_val<=test_val
    min_a=test_a;
    test_a=min_a/2+max_a/2;

    if abs(test_a-max(a))<d
        test_a=max(a);
        cond=false;
    end


    taus(ind)=test_a;
    [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
    test_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);
    else
        cond=false;
    end

end


end

% 
% cond=true;
% max_a=min(test_a+d,max(a));
% min_a=max(test_a-d,min(a));
% test_a=min_a/2+max_a/2;
% 
% while cond
% d=d/2;
% taus(ind)=min_a;
% [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
% left_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);
% 
% if left_val<=test_val
%     max_a=test_a;
%     test_a=min_a/2+max_a/2;
% 
%     if abs(test_a-min(a))<d
%         test_a=min(a);
%         cond=false;
%     end
% 
%     taus(ind)=test_a;
%     [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
%     test_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);
% 
% else 
% taus(ind)=max_a;
% [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
% right_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);
%     if right_val<=test_val
%     min_a=test_a;
%     test_a=min_a/2+max_a/2;
% 
%     if abs(test_a-max(a))<d
%         test_a=max(a);
%         cond=false;
%     end
% 
% 
%     taus(ind)=test_a;
%     [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
%     test_val=r(ind)+(P(ind,:)*v-v(ind))/h(ind);
%     else
%         cond=false;
%     end
% 
% end
% 
% if d<0.00001
%     cond=false;
% end
% 
% end


        mv=test_val;
        new_tau=test_a;
end

