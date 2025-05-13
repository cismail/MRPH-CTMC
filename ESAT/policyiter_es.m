function [tauVec]=policyiter_es(Q,mu,lagr,ind0,v,TauMat)

N=length(Q);
tauMax=8;
dt=0.05;
tauTVec=0:dt:tauMax;
L=length(tauTVec);
for ind1=1:L
    for ind2=1:L
        tau_vec=[tauTVec(ind1),tauTVec(ind2)];
        NVec=1:N;
        NVec(ind0)=[];
        TauMat(ind0,NVec)=tau_vec;
        [meanAoII(ind1,ind2),MeanSamplingRate(ind1,ind2),r,h,P]=M3RPhaseType(Q,mu,TauMat,lagr);
        val(ind1,ind2)=r(ind0)/h(ind0)+(P(ind0,:)*v-v(ind0))/h(ind0);

    end
end

[~,ix]=min(val(:));

[row,col] = ind2sub([L,L],ix);

tau_1=tauTVec(row);
tau_2=tauTVec(col);

tauVec=zeros(1,N);
tauVec(NVec)=[tau_1,tau_2];

end