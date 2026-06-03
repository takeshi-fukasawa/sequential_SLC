function [out,other_vars]=phat_update_func(phat,theta,...
    middle_vars,...
    beta,u1_0,u2_0)

    phat1=phat(:,:,1);
    phat2=phat(:,:,2);

    if isempty(middle_vars)==1
        middle_vars=construct_vars_from_phat_func(phat,...
        beta,u1_0,u2_0);
    end

    H1=middle_vars{1};
    H2=middle_vars{2};
    z1=middle_vars{3};
    z2=middle_vars{4};

    phat1_updated=phat1;
    phat2_updated=phat2;
    
    phat1_updated(:,2) = normcdf(H1*theta+z1);
    phat1_updated(:,1) = 1 - phat1_updated(:,2);
    phat2_updated(:,2) = normcdf(H2*theta+z2);
    phat2_updated(:,1) = 1 - phat2_updated(:,2);

    phat_updated=cat(3,phat1_updated,phat2_updated);

    %% Relaxation method
    weight_param=1;
    phat_updated=phat_updated*weight_param+phat*(1-weight_param);

    out={phat_updated};
    other_vars=[];

end
