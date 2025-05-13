function [age,cost,Table]=SMP_solver(Q,mu,target,opt)

sens=opt.lagr_sens;

lagr_max=opt.lagr_max;
lagr_min=opt.lagr_min;

old_lagr=-1;

if isequal(opt.lagr_solver,'bisection')
    cond=true;
    lagr_test=lagr_max/2+lagr_min/2;
    
    [obj,Table,age,cost] = SMP_lag_func_opt(Q,mu,lagr_test,opt);
    tic;
    while cond
        if cost<target
            lagr_max=lagr_test;
            lagr_test=lagr_max/2+lagr_min/2;
        else
            lagr_min=lagr_test;
            lagr_test=lagr_max/2+lagr_min/2;
        end
        [obj,Table,age,cost] = SMP_lag_func_opt(Q,mu,lagr_test,opt);
        cond=and(abs(cost-target)>sens,abs(old_lagr-lagr_test)>opt.lagr_sens);
        if lagr_test<sens
            lagr_test=0;
            [obj,Table,age,cost] = SMP_lag_func_opt(Q,mu,lagr_test,opt);
            cond=false;
        end
        if opt.display_iter
            disp(['Lagrangian:   ',num2str(lagr_test)])
            disp(Table)
        end
        t_end=toc;
        opt.tocs=[opt.tocs,t_end];
        old_lagr=lagr_test;

    end

elseif isequal(opt.lagr_solver,'search')

    ld=opt.lagr_stepsize;

    lagr_vec=lagr_min:ld:lagr_max;

    for i=1:length(lagr_vec)
        lagr_test=lagr_vec(i);
        [obj,Table,age,cost] = SMP_lag_func_opt(Q,mu,lagr_test,opt);
        cond=abs(cost-target)>sens;
        if not(cond) break; end
    end

else
    error('Undefined Lagrangian finder')
end


if opt.display_endtable
    disp(Table)
    disp(opt.tocs)
end


end