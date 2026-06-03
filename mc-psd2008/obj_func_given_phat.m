function [obj_val,grad]=obj_func_given_phat(theta,phat,...
middle_vars,...
beta,u1_0,u2_0,Data)

    if isempty(middle_vars)==1 
        middle_vars=construct_vars_from_phat_func(phat,...
        beta,u1_0,u2_0);
    end

     Data.H1 = middle_vars{1};
     Data.H2 = middle_vars{2};
     Data.z1 = middle_vars{3};
     Data.z2 = middle_vars{4};

    [obj_val,grad]=LL_NPL(theta,Data);

end
