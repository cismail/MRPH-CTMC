clear
addpath("SMPfunctions\")

% Q=[-1.2,0.7,0.3,0.1,0.1;
%     0.2,-0.8,0.4,0.1,0.1;
%     0.1,0.7,-1,0.1,0.1;
%     0.2,0.2,0.2,-0.8,0.2;
%     0.3,0.3,0.25,0.25,-1.1]/2;

mu=1;


Ls=[0.6,0.75]; %Q3

mu=1;
Q=[-Ls(1),Ls(1);Ls(2),-Ls(2)];

%
% N=5; %Q1
% sigma=0.7;
% Q=ones(N)/(N-1)-(1+1/(N-1))*eye(N);
% Q=sigma*Q;
%
%
% Q=[-102.5,100,2.5;
%      0.5,-7.5,7;
%      4,0.1,-4.1];
% Q=Q./100;
%

[Q]=QNN(.1,1,0.8,1.2,3);
N=length(Q);
mu=1;





opt.eps_value_dist=1e-2;
opt.eps_pol_dist=(N-1)*1e-2;
opt.lagr_max=1000;
opt.lagr_min=0;
opt.min_a=0;
opt.max_a=30; %30
opt.stepsize=1e-1;
opt.lagr_solver='bisection';
opt.policy_solver='search';
opt.initial_pol=opt.min_a*ones(1,N);
opt.display_iter=true;
opt.display_endtable=true;
opt.lagr_sens=1e-3;
opt.display_time=true;
opt.tocs=[];
opt.lagr_stepsize=1;
opt.stepsize_hyb=1;

[~,~,~,~,age_alw,cost_alw] = reward_lambda_tau(Q,mu,zeros(1,N),0);
[rewardN,~,cycleVecN,pikn,meanAoIIn,cost_nvr] = reward_lambda_tau(Q,mu,100*ones(1,N),0);

age_nvr=min(rewardN./cycleVecN);
Lims=[age_alw,age_nvr;cost_alw,cost_nvr];
disp(Lims)



rate=0.05:0.05:0.25;
for i=1:length(rate)
    target=rate(i);
    counter_bs=tic;
    opt.policy_solver='bisection';
    [age_bs(i),cost_bs(i),Table]=SMP_solver(Q,mu,target,opt);
    tau_bs=Table.Taus(end,:);
    time_bs(i)=toc(counter_bs);


    counter_st=tic;
    [tau_st(i),age_st(i),cost_st(i)] = ST_bs(Q,mu,target,50);
    t_st(i)=toc(counter_st);

    counter_pt=tic;
    [tau_pt(i),age_pt(i),cost_pt(i)] = PT_bs(Q,mu,target,50);
    t_pt(i)=toc(counter_pt);

%     counter_eat=tic;
%     [age_eat(i),cost_eat(i),Table]=SMP_solver(Q,mu,target,opt);
%     t_eat(i)=toc(counter_eat);
i
save REV_5
end

% 
% figure(51);
% subplot(2,1,1)
% hold on;
% grid on;
% plot(rate,age_bs,'b--','LineWidth',2)
% plot(rate,age_gdbs,'k--','LineWidth',2)
% xlabel('$\rho$','Interpreter','latex')
% ylabel('average AoII','Interpreter','latex')
% 
% subplot(2,1,2)
% hold on;
% grid on;
% plot(rate,time_bs,'b--','LineWidth',2)
% plot(rate,time_gdbs,'k--','LineWidth',2)
% xlabel('$\rho$','Interpreter','latex')
% ylabel('time','Interpreter','latex')
% 
% 
% figure(49);
% subplot(2,1,1)
% hold on;
% grid on;
% plot(rate,age_bs,'b--','LineWidth',2)
% plot(rate,age_sr,'r--','LineWidth',2)
% plot(rate,age_gd,'k--','LineWidth',2)
% xlabel('$\rho$','Interpreter','latex')
% ylabel('average AoII','Interpreter','latex')
% 
% subplot(2,1,2)
% hold on;
% grid on;
% plot(rate,time_bs,'b--','LineWidth',2)
% plot(rate,time_sr,'r--','LineWidth',2)
% plot(rate,time_gd,'k--','LineWidth',2)
% xlabel('$\rho$','Interpreter','latex')
% ylabel('time','Interpreter','latex')



% legend('bisection (anl.)','exhaustive search (anl.)','bisection (sim.)','exhaustive search (sim.)','always sampling bound','never sampling bound',fontsize=14)


% legend('$Q_1$','$Q_2$','$Q_3$','simulation','Interpreter','latex',fontsize=14)




%
% figure;
% hold on;
% grid on;
% plot(rate,age_bs,'ro','LineWidth',2)
% % plot(rate,age_bs,'bx','LineWidth',2)
% plot(rate,Perbs,'r--','LineWidth',2)
% % plot(rate,Persr,'b:','LineWidth',2)
%
% % plot(rate,Lims(1,1)*ones(size(rate)),'m--','LineWidth',2)
% % plot(rate,Lims(1,2)*ones(size(rate)),'k--','LineWidth',2)
% xlabel('sampling rate constraint')
% ylabel('average AoII')
% % legend('bisection (anl.)','exhaustive search (anl.)','bisection (sim.)','exhaustive search (sim.)','always sampling bound','never sampling bound',fontsize=14)
%
% [Aopt,Ropt,opt_tau] = CTMC_anl_hom_opt(N,sigma,mu,target);


figure;
grid on;
hold on;
plot(rate,age_pt,'k-.','LineWidth',2)
plot(rate,age_st,'r--','LineWidth',2)
plot(rate,age_bs,'b-','LineWidth',2)
plot(rate,ageESAT,'g:.','LineWidth',2)
xlim([rate(1) rate(end)])
legend('PS','ST','EAT',"ESAT","Interpreter", 'latex')
xlabel('$\rho$','Interpreter','latex')
ylabel('MAoII','Interpreter','latex')

