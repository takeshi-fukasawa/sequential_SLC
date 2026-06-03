function [indobs,indobs_a_1]=indobs_func(aobs_1,zexo_obs,z_exoval)

    [nobs,nplayer]=size(aobs_1);
    [numz_exo,dim_zexo]=size(z_exoval);
    
    dummy_zexo_obs = prod(reshape(zexo_obs,nobs,1,dim_zexo)==reshape(z_exoval,1,numz_exo,dim_zexo),3);%nobs*numz_exo
    indzexo = dummy_zexo_obs*[ 1:numz_exo ]';%nobs*1
    
    twop = kron(ones(nobs, 1), 2.^(nplayer-[1:nplayer]));
    indobs_a_1 = sum(aobs_1.*twop, 2)+1;
    indobs = (indzexo-1).*(2^nplayer) + indobs_a_1; % 1 to numz*(2^nplayer), for the states of the game
    
    end 
    