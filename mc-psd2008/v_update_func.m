function [out,other_vars]=v_update_func(v,theta,...
    middle_vars,...
    beta,u1_0,u2_0)

    v1=v(:,:,1);
    v2=v(:,:,2);

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

    v1_0_updated=u1_0+beta.*EV1_0;
    v1_1_updated=h1*theta+beta.*EV1_1;

    v2_0_updated=u2_0+beta.*EV2_0;
    v2_1_updated=h2*theta+beta.*EV2_1;

    v_updated=cat(3,[v1_0_updated,v1_1_updated],...
     [v2_0_updated,v2_1_updated]);

   
    %% Relaxation method
    weight_param=1;
    v_updated=v_updated*weight_param+v*(1-weight_param);

    out={v_updated};
    other_vars=[];

end
