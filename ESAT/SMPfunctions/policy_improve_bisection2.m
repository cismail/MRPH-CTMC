function [new_tau,mv] = policy_improve_bisection2(Q,mu,taus,lagr,v,a,d,ind)

min_a=min(a);
max_a=max(a);
test_a=min_a/2+max_a/2;

taus(ind)=test_a;
[r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
test_val=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);

cond=true;

while cond
taus(ind)=test_a-d;
[r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
left_val=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);

if left_val<=test_val
    max_a=test_a;
    test_a=min_a/2+max_a/2;

    if abs(test_a-min(a))<d
        test_a=min(a);
        cond=false;
    end

    taus(ind)=test_a;
    [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
    test_val=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);

else 
taus(ind)=test_a+d;
[r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
right_val=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);
    if right_val<=test_val
    min_a=test_a;
    test_a=min_a/2+max_a/2;

    if abs(test_a-max(a))<d
        test_a=max(a);
        cond=false;
    end


    taus(ind)=test_a;
    [r,P,h,~,~,~] = reward_lambda_tau(Q,mu,taus,lagr);
    test_val=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);
    else
        cond=false;
    end

end


end



        mv=test_val;
        new_tau=test_a;
end

