function [new_tau,mv] = policy_improve_search(Q,mu,taus,lagr,v,a,ind)

        for k=1:length(a)
            taus(ind)=a(k);
            [r,P,h,pik,meanAoII,MeanSamplingRate] = reward_lambda_tau(Q,mu,taus,lagr);
            val(k)=r(ind)/h(ind)+(P(ind,:)*v-v(ind))/h(ind);
%             [A(k),D(k),C(k),M3(k),res(k),dres(k)] = f98(Q,mu,v,ind,a(k),lagr);
%             val2(k)=(A(k)+lagr*C(k)+M3(k))/D(k);
        end
        
        [mv,ix]=min(val);
        new_tau=a(ix);


%         [new_tau_gd,mv_gd] = policy_improve_gdesc(Q,mu,taus,lagr,v,a,ind);
end

