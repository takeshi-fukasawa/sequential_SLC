    %% Create other_input_Phi_func, other_input_obj_func, x_initial0_cell, Phi_func, obj_func, param_initial_all, param_min, param_max
    
    %%%%%%%%%%%%
    %%% Generate data and solve equilibrium
    run DGP.m

    run solve_equil_for_seq.m

    run DGP_instruments.m


    sd=1;
    sigma_const=sigma_const_true.*2.*sd.*rand(1,nstart);
    sigma_x=sigma_x_true(:).*2.*sd.*rand(3,nstart);
    sigma_p=sigma_p_true.*2.*sd.*rand(1,nstart);
    sigma_vec_est=[sigma_const;sigma_x;sigma_p];
    %%sigma_vec_est=sigma_vec_est*10;%%%

    %V_initial0=-log(S_0t_data.*weight);
    delta_initial0=log(S_jt_data)-log(S_0t_data)-rho_est.*log(S_jt_given_g_data);% Initial value of delta


    Pr0_initial0=ones(1,ns,1,T);

    V_initial_global=V_initial0;


    other_input_Phi_func={v_est,X_jt,weight,S_jt_data,S_0t_data,...
            weight_V,x_V,beta_C,tune_param,Newton_spec};
    other_input_obj_func={weight,v_est,X_jt,X,Z,PX,W,...
            S_jt_data,S_0t_data,...
            weight_V,x_V,beta_C,tune_param,Newton_spec};



    param_initial_all=[sigma_vec_est([2,3,5],:)];%n_param*nstart

    x_initial0_cell={V_initial0};
    Phi_func=@V_param_update_func;
    obj_func=@obj_func_given_V;

    n_param=size(param_initial_all,1);
    param_min=zeros(n_param,1)-10;
    param_max=10*ones(n_param,1);
