function [Q]=QNN(s1,s2,q1,q2,N)

S=linspace(s1,s2,N);

for i=1:N
    q=S(i)*linspace(q1,q2,N-1)/(N-1);
    Q(i,:)=circshift([0,q],i-1);
end

E=Q*ones(N,1);
Q(1:N+1:end)=-E;


end