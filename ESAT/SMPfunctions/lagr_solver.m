function [age,cost] = lagr_solver(Q,mu,opt)

cond=true;

opt.lagr_max=lagr_max;
opt.lagr_min=lagr_min;
lagr_test=lagr_max/2+lagr_min/2;

[obj,Table,age,cost] = SMP_lag_func(Q,mu,lagr_test,opt);

if isequal(opt.lagr_solver,'bisection')
    a_min=opt.min_a;
    a_max=opt.max_a;
    d=opt.stepsize;
    a=[a_min,a_max];



    while cond
        if cost<target
            lagr_max=lagr_test;
            lagr_test=lagr_max/2+lagr_min/2;
        else
            lagr_min=lagr_test;
            lagr_test=lagr_max/2+lagr_min/2;
        end
        %disp(['Lagrangian:   ',num2str(lagr_test)])
        [obj,Table,age,cost] = SMP_lag_func(Q,mu,lagr_test,opt);
        cond=abs(cost-target)>sens;
    end

elseif isequal(opt.lagr_solver,'search')


elseif isequal(opt.lagr_solver,'gradient')

else
    error('Undefined policy solver')
end






end

