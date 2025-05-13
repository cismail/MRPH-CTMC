function [meanAoII,MeanSamplingRate,reward,h,T_ord]=M3RPhaseType(Q,mu,TauMat,lagr)


N=length(Q);

md_an=zeros(1,N);
ma_an=zeros(1,N);
mr_an=zeros(1,N);
Tan=zeros(1,N);
T_ord=zeros(N,N);

for i=1:3
Nvec=1:N;
Nvec(i)=[];
tr_order=[i,Nvec];


tau_i_vec=TauMat(i,:);
tau_i_vec(i)=[];

[tau_low,min_t]=min(tau_i_vec);
[tau_high,max_t]=max(tau_i_vec);

if tau_low==tau_high
    max_t=2;
end

Qi=Q(tr_order,tr_order);

Sigma=-diag(Qi);
i0=1;
Qt=Qi;
Qt(i0,:)=[];
Qt(:,i0)=[];
B0=Qt;
B1=Qt;
B1(min_t,min_t)=B1(min_t,min_t)-mu;

B2=B1;
B2(max_t,max_t)=B2(max_t,max_t)-mu;


Beta_1=Qi(i0,:)./Sigma(i0);
Beta_1(i0)=[];
Beta_2=Beta_1*expm(B0*tau_low);
Beta_3=Beta_2*expm(B1*(tau_high-tau_low));





f0=zeros(2,1); f_low=f0; 

f_low(min_t)=1;
f_high=ones(2,1);












H=[Sigma(1:2)+mu];





md_an(i)=Int_tx(Beta_1,B0,0,tau_low,[Qi(2,1);Qi(3,1)])+Int_tx(Beta_2,B1,tau_low,tau_high,[Qi(2,1);Qi(3,1)]+mu*f_low)+Int_txinf(Beta_3,B2,tau_high,[Qi(2,1);Qi(3,1)]+mu*f_high)+1/Sigma(1);
ma_an(i)=Int_tx22(Beta_1,B0,0,tau_low,[Qi(2,1);Qi(3,1)])+Int_tx22(Beta_2,B1,tau_low,tau_high,[Qi(2,1);Qi(3,1)]+mu*f_low)+Int_tx22inf(Beta_3,B2,tau_high,[Qi(2,1);Qi(3,1)]+mu*f_high);



Tan(1)=Int_1(Beta_1,B0,tau_low,[Qi(2,1);Qi(3,1)])+Int_1(Beta_2,B1,tau_high-tau_low,[Qi(2,1);Qi(3,1)])+Int_1inf(Beta_3,B2,[Qi(2,1);Qi(3,1)]);
Tan(1+min_t)=Int_1(Beta_2,B1,-tau_low+tau_high,f_low)+Int_1inf(Beta_3,B2,f_low);
Tan(1+max_t)=Int_1inf(Beta_3,B2,(f_high-f_low));







if min_t==1

    mr_low2an=Int_1(Beta_2,B1,tau_high-tau_low,[0;1])*Qi(max_t+1,min_t+1)+Beta_2*[1;0];
    mr_high2an=Int_1inf(Beta_3,B2,[1;0])*(Qi(min_t+1,max_t+1))+Int_1inf(Beta_3,B2,[0;1])*Qi(max_t+1,min_t+1)+Beta_3*[0;1];

else

    mr_low2an=Int_1(Beta_2,B1,tau_high-tau_low,[1;0])*Qi(max_t+1,min_t+1)+Beta_2*[0;1];
    mr_high2an=Int_1inf(Beta_3,B2,[0;1])*(Qi(min_t+1,max_t+1))+Int_1inf(Beta_3,B2,[1;0])*Qi(max_t+1,min_t+1)+Beta_3*[1;0];
   
end

% MAoII(i)=ma_an/md_an;
% 
% R(i)=(mr_low2an+mr_high2an)/double(md_an);

mr_an(i)=(mr_low2an+mr_high2an);

T_ord(i,[i,Nvec])=Tan;

end

e=ones(N,1);

% steady-state solution of the embedded chain
pik=e'*inv(T_ord+e*e'-eye(N,N));




num=0;
den=0;
for i=1:N
    num = num+pik(i)*ma_an(i);
    den=den+pik(i)*md_an(i);
end

meanAoII = num/den;


sum1=0;
for i=1:N
    sum1 =sum1 + pik(i)*mr_an(i);
end

MeanSamplingRate=sum1/den;

reward=ma_an+lagr*mr_an;
h=md_an;


end