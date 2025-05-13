function [reward,P,cycleVec,pik,meanAoII,MeanSamplingRate] = reward_lambda_tau(Q,mu,tauVec,lambda)


N=length(Q);
if length(tauVec)==1
    tauVec=tauVec*ones(1,N);
end


sigmaVec=-diag(Q);
Qt=Q+diag(sigmaVec);
QtEmbedded=Qt;

for k=1:N
    QtEmbedded(k,:)=QtEmbedded(k,:)/sum(QtEmbedded(k,:));
    p1=QtEmbedded(k,:);
    c1 = cumsum([0,p1(:).']);
    c1 = c1/c1(end); 
    QtEmbeddedCum(k,:) = c1;
end



e=ones(N,1);

eN=ones(N,1);
e2N=ones(2*N-1,1);
meanVec=zeros(N,1);
areaVec=zeros(N,1);
sampleVec=zeros(N,1);
% now we are going to establish the discrete-time MC embedded at the epochs
% of synchronization

P=zeros(N,N);




 for i=1:N
tau=tauVec(i);

    beta=Q(i,:)./sigmaVec(i);
    beta(i)=0;
    beta(N+1)=0;
    
    A=[Q,Q(:,i)];
    A(i,N+1)=0;
    A(:,i)=zeros(N,1);
    A(i,i)=Q(i,i);
    A(N+1,:)=zeros(1,N+1);
panl=beta*expm(A*tau);
B=A(1:N,1:N);
b=A(1:N,N+1);


e=ones(N,1);

po=panl(N+1);
alpha=beta(1:N);

mid_term_d=-expm(B*tau)*tau+expm(B*tau)*inv(B)-inv(B);
DurT=alpha*mid_term_d*e/po;


mid_term_a=expm(B*tau)*  (tau^2*B^(-1)-2*B^(-2)*tau+2*B^(-3))  -2*B^(-3);
AreaT=alpha*mid_term_a*b/po/2;

if isnan(DurT) DurT=0; end
if isnan(AreaT) AreaT=0; end


eN=ones(2*height(Q),1);
Cbig=[Qt,eye(N)*mu];
Cbig(:,N+i)=Qt(:,i);
Cbig(i,N+i)=0;
Cbig(:,i)=zeros(N,1);
Cbig(N+1:2*N,:)=zeros(N,2*N);
Cbig=Cbig - diag(Cbig*eN);
C=Cbig(1:N,1:N);
c=mu*ones(N,1);
c(i)=0;

alpc=panl(1:N);

d_a=-alpc*inv(C)*e;
d=d_a+DurT*po+(1-po)*tau+1/sigmaVec(i);

a_a=alpc*C^(-2)*e;

a=a_a+(1-po)*tau^2/2+tau*d_a+po*AreaT;
cycleVec(i)=d;
areaVec(i)=a;

    trans_rate=mu./(sigmaVec+mu);


  E(1:N,1:N)=Qt./sigmaVec.*(1-trans_rate);
  E(i,:)=zeros(N,1);

    TravelCost=inv(eye(N)-E);
    TravelCost(:,i)=zeros(N,1);

    
    sampleVec(i)=alpc*TravelCost*ones(N,1); 
    
    alpha=zeros(1,2*N-1);
    alpha(i)=1;
    for j=1:N
        P(i,j)=-alpc*inv(C)*Cbig(1:N,N+j);
    end
    P(i,i)=P(i,i)+po;

    reward(i)=a+lambda*sampleVec(i);

 end




% steady-state solution of the embedded chain
pik=e'*inv(P+e*e'-eye(N,N));
if any(pik<0)
    warning('non-positive pik')
end
pik=max(pik,zeros(size(pik)));


num=0;
den=0;
for i=1:N
    num = num+pik(i)*areaVec(i);
    den=den+pik(i)*cycleVec(i);
end

meanAoII = num/den;


sum1=0;
for i=1:N
    sum1 =sum1 + pik(i)*sampleVec(i);
end

MeanSamplingRate=sum1/den;
end

