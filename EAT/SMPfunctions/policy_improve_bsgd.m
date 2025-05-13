function [new_tau,mv] = policy_improve_bsgd(Q,mu,lagr,v,a,d,ind)

max_a=max(a);

N=length(Q);
e=ones(N-1,1);
sigma=-Q(ind,ind);

Qi=Q;
Qi(ind,:)=[];
Qi(:,ind)=[];

beta=-Q(ind,:)/Q(ind,ind);
beta(ind)=[];

A1=Qi;
A2=Qi-mu*eye(N-1);

D=zeros(N-1);
for i=1:N-1
    for j=1:N-1
        if i~=j
            D(i,j)=-A2(i,j)/A2(i,i);
        end
    end
end
F=inv(eye(N-1)-D);

cd=-beta*inv(A1)*e+1/sigma;
ca=beta*(A1)^(-2)*e;

M1=inv(A1)-inv(A2);
M2=A1^(-2)-A2^(-2);

tau=sym('tau');

betax=beta*expm(tau*A1);
B=eye(N-1)*mu;
% P=-betax*inv(A2)*B;

vi=v;
vi(ind)=[];
V=diag(vi-v(ind));

M3=-inv(A2)*B*V;

% pv=P*V*e;
% a=betax*(tau*M1-M2)*e+ca;
% d=beta*expm(tau*A1)*M1*e+cd;
% c=betax*F*e;
% res1=(a+lagr*c+pv)/d;

res_n=betax*(tau*M1-M2+lagr*F+M3)*e+ca;
res_d=betax*M1*e+cd;
res=res_n/res_d;

f=res_n;
df=betax*( M1*A1*tau+M1-A1*M2+lagr*A1*F+ A1*M3)*e;
g=res_d;
dg=betax*A1*M1*e;

dres=(g*df-dg*f)/g^2;


% eta=1;
% old_tau=0;
cond=true;

min_a=min(a);
max_a=max(a);
test_a=min_a/2+max_a/2;


        while cond
        
                    dt=subs(dres,test_a);

                    if dt>0
                        max_a=test_a;
                        test_a=min_a/2+max_a/2;
                    else
                        min_a=test_a;
                        test_a=min_a/2+max_a/2;
                    end






    if abs(test_a-min(a))<d
        test_a=min(a);
        cond=false;
    elseif abs(test_a-max(a))<d
        test_a=max(a);
        cond=false;
    elseif abs(max_a-min_a)<d
        cond=false;
    end                    


                    
        
        
        end


mv=double(subs(res,test_a));
new_tau=test_a;
end

