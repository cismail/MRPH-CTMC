function [new_tau,mv] = policy_solver(Q,mu,lagr,taus,v,ind,opt)

if isequal(opt.policy_solver,'bisection')
    a_min=opt.min_a;
    a_max=opt.max_a;
    d=opt.stepsize;
    a=[a_min,a_max];

    [new_tau,mv] = policy_improve_bisection2(Q,mu,taus,lagr,v,a,d,ind);

elseif isequal(opt.policy_solver,'search')
    a_min=opt.min_a;
    a_max=opt.max_a;
    d=opt.stepsize;
    a=a_min:d:a_max;
    [new_tau,mv] = policy_improve_search(Q,mu,taus,lagr,v,a,ind);

elseif isequal(opt.policy_solver,'gradient')


elseif isequal(opt.policy_solver,'hybrid')
    a_min=opt.min_a;
    a_max=opt.max_a;
    d=opt.stepsize_hyb;
    a=a_min:d:a_max;
    [new_tau,mv] = policy_improve_search(Q,mu,taus,lagr,v,a,ind);
    if new_tau>d
    hyb_a=[new_tau-d,new_tau+d];
    d_bs=opt.stepsize;
    [new_tau,mv] = policy_improve_bisection2(Q,mu,taus,lagr,v,a,d_bs,ind);
    end

else
    error('Undefined policy solver')
end


end

