function [obj_val,grad]=obj_func_given_v(theta,v,...
middle_vars,...
beta,u1_0,u2_0,Data)

    if isempty(middle_vars)==1 
        middle_vars=construct_vars_from_v_func(v,...
        beta,u1_0,u2_0);
    end

    EV1_0=middle_vars{1};
    EV1_1=middle_vars{2};
    EV2_0=middle_vars{3};
    EV2_1=middle_vars{4};
    h1=middle_vars{5};
    h2=middle_vars{6};

     H1 = h1;
     z1 = -u1_0+beta.*(EV1_1-EV1_0);

     H2 = h2;
     z2 = -u2_0+beta.*(EV2_1-EV2_0);

     Data.H1 = H1;
     Data.H2 = H2;
     Data.z1 = z1;
     Data.z2 = z2;

    [obj_val,grad]=LL_NPL(theta,Data);

end
