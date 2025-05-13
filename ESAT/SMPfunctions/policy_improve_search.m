function [new_tau,mv] = policy_improve_search(Q,mu,taus,lagr,v,a,ind)

        for k=1:length(a)
            taus(ind)=a(k);
            [r,P,h,pik,meanAoII,MeanSamplingRate] = reward_lambda_tau(Q,mu,taus,lagr);
            val(k)=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);
        end
        
        [mv,ix]=min(val);
        new_tau=a(ix);
end

