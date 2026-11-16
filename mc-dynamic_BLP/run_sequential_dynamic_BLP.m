
file_name_create_variables="create_variables_dynamic_BLP.m";


n_param=3;
n_sim=n_market;
%%%%%%%
file_name_save_variables=append(save_path,"dynamic_BLP_seq_results");

if run_nfxp_spec==0
   file_name_save_variables=append(file_name_save_variables,"_SLC_only");
end
if spectral_spec==1
    file_name_save_variables=append(file_name_save_variables,"_spectral");
end
if SLC_tol==1e-4
    file_name_save_variables=append(file_name_save_variables,"_SLC_loose_tol");
end
if nstart>=10
    file_name_save_variables=append(file_name_save_variables,"_many_starting_values");
end

if bad_init_val_spec==1
    file_name_save_variables=append(file_name_save_variables,"_bad_init");
end
if SLC_line_search_spec==1
    file_name_save_variables=append(file_name_save_variables,"_line_search");
end

file_name_save_variables=append(file_name_save_variables,".mat");

SLC_positive_param_update_spec=0;
initial_x_solve_spec=1;%If initial_x_solve_spec==1, solve for the initial x (e.g., value function) that is consistent with the initial parameter guess. 
%%If initial_x_solve_spec==0, use the initial x specified in x_initial0_cell without solving for consistency.

%%check_conv_to_same_param_spec=0;

if check_conv_to_same_param_spec==1
    run_nfxp_spec=0;
    nstart=1;
    load(append(save_path,"dynamic_BLP_seq_results.mat"),...
    'param_nfxp_all','obj_val_nfxp_all')
     
    [min_val, temp] = min(obj_val_nfxp_all, [], 2);
    param_nfxp_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_nfxp_all_i = param_nfxp_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_nfxp_all_i), rows, temp);
        param_nfxp_est(:,i) = param_nfxp_all_i(idx);
    end

end

market_id=1;

    obj_val_SLC_all=NaN(n_sim,nstart);
    times_SLC_all=NaN(n_sim,nstart);
    n_iter_SLC_all=NaN(n_sim,nstart);
    count_G_eval_SLC_all=NaN(n_sim,nstart);
    count_obj_eval_SLC_all=NaN(n_sim,nstart);

    obj_val_nfxp_all=NaN(n_sim,nstart);
    times_nfxp_all=NaN(n_sim,nstart);
    n_iter_nfxp_all=NaN(n_sim,nstart);
    count_G_eval_nfxp_all=NaN(n_sim,nstart);
    count_obj_eval_nfxp_all=NaN(n_sim,nstart);
    
    param_SLC_all=NaN(n_market,nstart,n_param);
    param_nfxp_all=NaN(n_market,nstart,n_param);

    iter_info_SLC_all=[];
    results_nfxp_all=[];

for draw=1:n_sim%1:n_sim
    rng(100+draw);
    

    %% Create other_input_Phi_func, other_input_obj_func, x_initial0_cell, Phi_func, obj_func, param_initial_all, param_min, param_max

    %%% Generate data and solve equilibrium
    run DGP.m

    run solve_equil_for_seq.m

    run DGP_instruments.m
    
    mu_ijt_est=mu_func(sigma_vec_est,v_est,X_jt);
    mu_ijt_est=reshape(mu_ijt_est,J,ns,G,T)./eps_sd;


    V_initial0=zeros(1,ns,1,n_dim_V);%1*ns*1*n_dim_V
    V_initial_global=V_initial0;


    other_input_Phi_func={v_est,X_jt,weight,S_jt_data,S_0t_data,...
            weight_V,x_V,beta_C,tune_param,Newton_spec};
    other_input_obj_func={weight,v_est,X_jt,X,Z,PX,W,...
            S_jt_data,S_0t_data,...
            weight_V,x_V,beta_C,tune_param,Newton_spec};


    param_true=[sigma_vec_true(2);sigma_vec_true(3);sigma_vec_true(5)];

    x_initial0_cell={V_initial0};
    Phi_func=@V_param_update_func;
    obj_func=@obj_func_given_V;


    param_min=param_true*0-5;
    param_max=param_true*0+5;
    param_max(3)=0.5;%%% Price: Large variation; Setting too large sigma_p => Numerically unstable...

    rho_val=0;

    opt_options = [];

    middle_vars_spec_cell=[];

    spec_sequential=[];
    spec_sequential.run_nfxp=run_nfxp_spec;
    spec_sequential.spectral_spec=spectral_spec;
    spec_sequential.x_min_cell={param_min,[],[],[],[]};
    spec_sequential.x_max_cell={param_max,[],[],[],[]};
    spec_sequential.SLC_positive_param_update_spec=SLC_positive_param_update_spec;
    spec_sequential.initial_x_solve_spec=initial_x_solve_spec;
    spec_sequential.run_EPL_spec=run_EPL_spec;
    spec_sequential.line_search_spec=SLC_line_search_spec;

    spec_sequential.TOL=SLC_tol;

    n_vars=size(x_initial0_cell,2);

    spec_gmres=[];
    spec_gmres.ITER_MAX=400;
    spec_gmres.TOL=1e-8;

    %% Random initial parameters
    sd=1;
    if bad_init_val_spec==0
        sigma_const=sigma_const_true.*2.*sd.*rand(1,nstart);
        sigma_x=sigma_x_true(:).*2.*sd.*rand(3,nstart);
        sigma_p=sigma_p_true.*1.5.*sd.*rand(1,nstart);
    else%bad_init_val_spec==1
        sigma_const=sigma_const_true.*3.*sd.*ones(1,nstart);
        sigma_x=sigma_x_true(:).*3.*sd.*ones(3,nstart);
        sigma_p=sigma_p_true.*1.5*sd.*ones(1,nstart);
    end
    
    sigma_vec_est=[sigma_const;sigma_x;sigma_p];%5*nstart

    for r=1:nstart
    %% Random draw of initial parameters
        param_initial=[sigma_vec_est(2,r);sigma_vec_est(3,r);sigma_vec_est(5,r)];
        
        if check_conv_to_same_param_spec==1
            %% Perturb NFXP sol and check SLC convergence
            %param_initial=[0;param_nfxp_est(draw,1);param_nfxp_est(draw,2);0;...
            %    param_nfxp_est(draw,3)]*1.001;
            param_initial=reshape(param_nfxp_est(draw,:),n_param,1)*1.001;
        end

        [param_SLC_all(draw,r,:),iter_info_SLC,X_sol_SLC,...
        iter_info_x,other_vars,speed_test_info,results_nfxp]=...
        run_sequential_est_func0(param_initial,x_initial0_cell{:},Phi_func,obj_func,middle_vars_spec_cell,...
        other_input_Phi_func,other_input_obj_func,...
        rho_val,opt_options,spec_sequential,spec_gmres,n_vars);

        obj_val_SLC_all(draw,r)=iter_info_SLC.obj_val;
        iter_SLC_all(draw,r)=iter_info_SLC.n_iter;
        times_SLC_all(draw,r)=iter_info_SLC.t_cpu;
        count_G_eval_SLC_all(draw,r)=iter_info_SLC.count_G_eval;
        count_obj_eval_SLC_all(draw,r)=iter_info_SLC.count_obj_eval;
        iter_info_SLC_all{draw,r}=iter_info_SLC;

        if isempty(results_nfxp)==0
            param_nfxp_all(draw,r,:)=results_nfxp.param_nfxp;

            iter_info_nfxp=results_nfxp.iter_info_nfxp;
            obj_val_nfxp_all(draw,r)=iter_info_nfxp.fval_NFXP;
            iter_nfxp_all(draw,r)=iter_info_nfxp.output_opt_NFXP.iterations;
            times_nfxp_all(draw,r)=iter_info_nfxp.t_nfxp;
            count_G_eval_nfxp_all(draw,r)=iter_info_nfxp.count_Phi_eval;
            count_obj_eval_nfxp_all(draw,r)=iter_info_nfxp.count_obj_eval;
            results_nfxp_all{draw,r}=results_nfxp;
        end
    end% iteration wrt r
    if save_spec==1 & check_conv_to_same_param_spec==0
        save(file_name_save_variables);
    end
end%draw

