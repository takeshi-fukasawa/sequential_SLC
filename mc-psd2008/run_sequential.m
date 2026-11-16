clear all

for draw=1:5
    rng(draw);
    
    
    nsims = 1;
    
    eqm_dgp=2;
    
    nobs=250;
    % Generate data from the model
    generate_data
    
    %% EPL iteration using generated data
    EPL_not_update_param_spec=0;
    initial_est_spec=1;% Use NPL first stage estimates
    % Estimate the model with consistent CCPs
    c = 1.0;
    nstart = 1;
    psd2008_estimate
    
    v_EPL_sol=v;
    param_EPL_sol=theta_EPL;
    
    v0=v;
    
    param_initial=theta0*1.0;%%%%
    param_initial=theta_NPL_1step_all(:,r,j)*1;%%%%
    v_init=v_EPL_sol*1.1;
    param_initial=param_EPL_sol*1.1;
    param_initial=theta_NPL_1step_all(:,r,j).*10;
    %param_initial=zeros(size(theta_NPL_1step_all(:,r,j)));
    v_init=5.*ones(size(v_EPL_sol));
    
    
    x_initial=reshape(v_init,4,2,2);
    
    Phi_func=@v_update_func;
    obj_func=@obj_func_given_v;
    middle_vars_func=@construct_vars_from_v_func;
    
    other_input_Phi_func={[],beta,u1_0,u2_0};
    other_input_obj_func={[],beta,u1_0,u2_0,Data_sim};
    
    other_input_middle_vars_func={beta,u1_0,u2_0};
    
    middle_vars_spec_cell.middle_vars_func=middle_vars_func;
    middle_vars_spec_cell.other_input_middle_vars_func=other_input_middle_vars_func;
    
    opt_options = [];
    
    spec_sequential=[];
    spec_sequential.spectral_spec=0;
    spec_sequential.line_search_spec=0;
    
    n_vars_x=1;
    
    spec_gmres=[];
    spec_gmres.ITER_MAX=100;
    spec_gmres.TOL=1e-8;
    rho_val=0;
    
    [param_SLC,iter_info_SLC,X_sol_SLC,...
            iter_info_x,other_vars,speed_test_info,results_nfxp,results_EPL]=...
        run_sequential_est_func0(param_initial*1.0,x_initial,Phi_func,obj_func,middle_vars_spec_cell,...
    other_input_Phi_func,other_input_obj_func,...
    rho_val,opt_options,spec_sequential,spec_gmres,n_vars_x);
    
    v_SLC=X_sol_SLC{2};
end


