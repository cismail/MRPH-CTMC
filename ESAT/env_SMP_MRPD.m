%clear
addpath("SMPfunctions\")
addpath("fs\")


mu=100;
Q=[-102.5,100,2.5; %Q2
     5,-75,70;
     40,1,-41];
Q=Q./100;

N=length(Q);







opt.eps_value_dist=1e-2;
opt.eps_pol_dist=(N-1)*1e-2;
opt.lagr_max=10;
opt.lagr_min=0;
opt.min_a=0;
opt.max_a=10; %30
opt.stepsize=1e-1;
opt.lagr_solver='bisection';
opt.policy_solver='bisection';
opt.initial_pol=TauMat0;
opt.display_iter=true;
opt.display_endtable=true;
opt.lagr_sens=1e-3;
opt.display_time=true;
opt.tocs=[];
opt.lagr_stepsize=1;
opt.stepsize_hyb=1;







rate=0.1:0.1:1;
for i=1:length(rate)

target=rate(i);

counter_esat=tic;
[ageESAT(i),costESAT(i),Table]=SMP_solver(Q,mu,target,opt);
t_esat(i)=toc(counter_esat);
% 
% tau_bs=Table.Taus(end,:);
% %[Perbs(i),Acostbs(i)] = env_func_taus(Q,tau_bs,mu,1e5,0);
% TauMat(:,:,i)=Table.Taus(end-N+1:end,:);
% [meanAoII_ESAT(i),MeanSamplingRate_ESAT(i)]=M3RPhaseType_sim(Q,mu,TauMat(:,:,i),MaxIND);
% save Q3_ESAT2
i
end

% 
% figure
% plot(rate,ageESAT,'b-')
% hold on
% plot(rate,meanAoII_ESAT,'bo')

% target=mean(Lims(2,:))
% rate=0.1:0.1:1;
% for i=1:length(rate)
% 
% target=rate(i);
% [age_bs(i),cost_bs(i),Table]=SMP_solver(Q,mu,target,opt);
% tau_bs=Table.Taus(end,:);
% [Perbs(i),Acostbs(i)] = env_func_taus(Q,tau_bs,mu,1e5,0);
% 
% end

