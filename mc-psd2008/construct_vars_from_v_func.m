function out=construct_vars_from_v_func(v,...
    beta,u1_0,u2_0)

    v1=v(:,:,1);
    v2=v(:,:,2);
    v1diff=v1(:,2)-v1(:,1);
    v2diff=v2(:,2)-v2(:,1);

    phat1=[1-normcdf(v1diff),normcdf(v1diff)];
    phat2=[1-normcdf(v2diff),normcdf(v2diff)];

    nstates=size(phat1,1);

            %u_i = h_i*theta
            h1 = [ones(nstates,1), phat2(:,2), [1;1;0;0]];
            h2 = [ones(nstates,1), phat1(:,2), [1;0;1;0]];
            
            %conditional transition matrices
            F1_0 = [phat2(:,1), phat2(:,2), zeros(nstates,2)]; %agent 1, a = 0
            F1_1 = [zeros(nstates,2), phat2(:,1), phat2(:,2)]; %agent 1, a = 1
            F2_0 = [phat1(:,1), zeros(nstates,1), phat1(:,2), zeros(nstates,1)]; %agent 2, a = 0
            F2_1 = [zeros(nstates,1), phat1(:,1), zeros(nstates,1), phat1(:,2)]; %agent 2, a = 1
            
   normsurplus_v1=normsurplus(v1);
   normsurplus_v2=normsurplus(v2);

    EV1_0=F1_0*normsurplus_v1;
    EV1_1=F1_1*normsurplus_v1;

    EV2_0=F2_0*normsurplus_v2;
    EV2_1=F2_1*normsurplus_v2;


       out={EV1_0,EV1_1,EV2_0,EV2_1,h1,h2};
end
