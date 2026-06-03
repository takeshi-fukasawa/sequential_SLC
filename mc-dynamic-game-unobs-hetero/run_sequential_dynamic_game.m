
global flag_vec relres_cell iter_gmres_cell resvec_cell

z_exoval=z_exo_val;
ptranz_exo_est=ptranz_exo;
theta=theta0;
v=v0;
bdisc=disfact_est;

%%%%%%%%%%%
v_spec=1;
%%%%%%%%%%%%%%


  run preliminary_construct_variables.m

  nobs_temp = size(aobs, 1)/T;

spec_sequential=[];

if unobs_hetero_spec==1 &1==1
    theta0=theta;
    %%% Recompute Initial CCP when unobs heterogeneity exists
    %%%% Use CCP estimates assuming z_exo2 does not exists
    ignore_zexo2_spec=1;%%% Even when zexo2 is unobserved, treat the variable as if observed.
    [v,prob0SP,theta_est]=estimate_initial_CCP_func(aobs,aobs_1,zexo_obs,mstate,ignore_zexo2_spec,...
        nobs_temp*T,nplayer,n_state,numz_exo_vec);
    theta=theta_est;
    %theta=theta0;%%%%
end


%%%%%%%%
dummy_aobs=(reshape(aobs,T,nobs_temp,nplayer,1)==reshape([0,1],1,1,1,naction));%T*nobs_temp*nplayer*naction

indobs_state_obs=reshape(indobs_a_1,T,nobs_temp);

if size(z_exoval,2)==2
    dummy_temp=...
        (reshape(zexo_obs(:,1),nobs_temp*T,1)==reshape(z_exoval(1:numz_exo_vec(1),1),1,numz_exo_vec(1)));
    indobs_zexo_obs=dummy_temp*[1:numz_exo_vec(1)]';
    indobs_state_obs = (indobs_zexo_obs-1).*(2^nplayer) + indobs_a_1;
end

indobs_state_obs_temp=indobs_state_obs(:);
for i=1:T*nobs_temp
    dummy_xobs_obs(i,indobs_state_obs_temp(i))=1;
end

%dummy_xobs_obs=sparse(dummy_xobs_obs);



if unobs_hetero_spec==0
    param_initial=theta;
    if v_spec==0
       %%% Use policy iteration mapping (cf. AM2007)
        pchoice = zeros(n_state, nplayer);
    
        % Calculate initial choice probabilities
        for i = 1:nplayer
            vdiff = v(:,i,2) - v(:,i,1);
            expvdiff = exp(vdiff);
            pchoice(:,i) = expvdiff ./ (1 + expvdiff);
        end
    
        x_initial_cell={pchoice};
        Phi_func=@pchoice_update_func;
        obj_func=@obj_func_given_pchoice;
        middle_vars_func=@compute_e_u_from_pchoice_func;
    elseif v_spec==1
       %%% Use the mapping on choice-specific value function
        x_initial_cell={v};
        Phi_func=@v_update_func;
        obj_func=@obj_func_given_v;
        middle_vars_func=@compute_e_u_from_v_func;
    end
    
    other_input_Phi_func={[],...
    ptranz_exo,bdisc,aval,mstate,ptranz_exo_kron};
    
    other_input_obj_func={[],...
    ptranz_exo, bdisc,...
    aval,mstate,dummy_count,ptranz_exo_kron};

    n_vars=size(x_initial_cell,2);

elseif unobs_hetero_spec==1
    
    param_initial=[theta(:);0];
    
    psteady = (1/n_state) * ones(n_state, 1);
    psteady_1=psteady(1:end-1);

    x_initial_cell={v,psteady_1};
    Phi_func=@v_psteady_update_func;
    obj_func=@obj_func_given_v_psteady_unobs_hetero;
    
    n_vars=size(x_initial_cell,2);

    other_input_Phi_func={[],numz_exo_vec,...
        bdisc,aval,mstate};
        
    other_input_obj_func={[],...
        numz_exo_vec,bdisc,aval,mstate,dummy_aobs,dummy_xobs_obs,indobs_state_obs};
        
end%unobs_hetero_spec==0 or 1


x_initial0_cell=x_initial_cell;
    
opt_options=[];

middle_vars_spec_cell=[];


spec_gmres=[];
spec_gmres.ITER_MAX=500;
spec_gmres.TOL=1e-6;


theta_init=theta;


spec_sequential.ITER_MAX=500;
if run_nfxp_spec==1
    spec_sequential.run_nfxp=1;
else
    spec_sequential.run_nfxp=0;
end

%spec_sequential.update_spec=0;
%spec_sequential.line_search_spec=1;
spec_sequential.common_spectral_coef_spec=0;

spec_sequential.compute_z_y_func=@compute_z_y_func_given_v_psteady;
%spec_sequential.Newton_update_func=@compute_z_func_given_v_psteady_unobs_hetero;

spec_sequential.analytical_spec=analytical_spec;%1;

spec_sequential.spectral_spec=spectral_spec;

%spec_sequential.nfxp_Newton_spec=1;
spec_sequential.nfxp_bhhh_spec=0;

if spec_sequential.analytical_spec==1
    spec_sequential.nfxp_bhhh_spec=0;%%%%%%
end

spec_sequential.initial_x_solve_spec=0; % If initial_x_solve_spec==1, solve for the initial x (e.g., value function) that is consistent with the initial parameter guess. 
% If initial_x_solve_spec==0, use the initial x specified in x_initial0_cell without solving for consistency.

%%spec_sequential.Anderson_acceleration=1;
%spec_sequential.SQUAREM_spec=1;

n_param=size(param_initial(:),1);
param_min=[];
param_max=[];
if 1==1
    param_min=-10*ones(size(param_initial));
    param_max=10*ones(size(param_initial));
end

if unobs_hetero_spec==1 & 1==1
    param_min(end)=0.0;
    param_max(end)=1.00;
    param_min(nplayer+2)=0;% To uniquely identify the parameter
    param_min(end-2)=0;%% More competitors=> Smaller profit
end

if isempty(param_min)==0 && isempty(param_max)==0
    spec_sequential.x_min_cell=cell(1,2*n_vars+1);
    spec_sequential.x_max_cell=cell(1,2*n_vars+1);

    spec_sequential.x_min_cell{1}=param_min;
    spec_sequential.x_max_cell{1}=param_max;
end

    spec_sequential.spec_x_update=[];
    spec_sequential.spec_x_update.ITER_MAX=1;
    spec_sequential.qNPL_spec=0;
    
        
%rho_val=0.1;
rho_val=0;
%%save_spec=0;

param_initial=param_initial(:);
param_initial_all=param_initial;
if nstart_seq>=2
    param_initial_all=[param_initial_all,...
        param_initial.*(0.5+rand(n_param,nstart_seq-1))];
end

param_initial_all(end,:)=rand(1,nstart_seq);
param_initial_all(end-3,:)=2*rand(1,nstart_seq);

%if unobs_hetero_spec==1
%    param_initial_all=3*rand(n_param,nstart_seq).*[theta_true(:);pi_param_true];
%end

param_initial_all=max(min(param_initial_all,param_max(:)),param_min(:));


    for r=1:nstart_seq
        
        if draw==1 & r==1
            param_SLC_all=NaN(n_sim,nstart_seq,n_param);
            param_nfxp_all=NaN(n_sim,nstart_seq,n_param);
        end
        [param_SLC_all(draw,r,:),iter_info_SLC ,X_sol_SLC{r},...
        iter_info_x,other_vars,speed_test_info,results_nfxp]=...
        run_sequential_est_func0(param_initial_all(:,r)*1.0,x_initial0_cell{:},Phi_func,obj_func,middle_vars_spec_cell,...
        other_input_Phi_func,other_input_obj_func,...
        rho_val,opt_options,spec_sequential,spec_gmres,n_vars);

        obj_val_SLC_all(draw,r)=iter_info_SLC.obj_val;
        iter_SLC_all(draw,r)=iter_info_SLC.n_iter;
        times_SLC_all(draw,r)=iter_info_SLC.t_cpu;
        count_G_eval_SLC_all(draw,r)=iter_info_SLC.count_G_eval;
        count_obj_eval_SLC_all(draw,r)=iter_info_SLC.count_obj_eval;
        
        if isempty(results_nfxp)==0
            param_nfxp_all(draw,r,:)=results_nfxp.param_nfxp;

            obj_val_nfxp_all(draw,r)=results_nfxp.iter_info_nfxp.fval_NFXP;
            iter_nfxp_all(draw,r)=results_nfxp.iter_info_nfxp.output_opt_NFXP.iterations;
            times_nfxp_all(draw,r)=results_nfxp.iter_info_nfxp.t_nfxp;
            count_G_eval_nfxp_all(draw,r)=results_nfxp.iter_info_nfxp.count_Phi_eval;
            count_obj_eval_nfxp_all(draw,r)=results_nfxp.iter_info_nfxp.count_obj_eval;
        end


    end% iteration wrt r
    
    if save_spec==1
        save_path="C:/Users/fukas/Dropbox/sequential/output/";
        file_name=append(save_path,"dynamic_game_unobs_hetero_seq_results_");
        if analytical_spec==1 
            if VFI_spec==2
                file_name=append(file_name,"_analytical_diff_EVFI.mat");
            else
                file_name=append(file_name,string(nobs),"_analytical_diff.mat");
            end
        else
            if VFI_spec==2
                file_name=append(file_name,"numerical_diff_EVFI.mat");
            else%VFI_spec==0
                file_name=append(file_name,"numerical_diff.mat");
            end
        end
        save(file_name)
    end

    tic
    for test_id=1:100
        obj_func(param_initial,x_initial_cell{:},other_input_obj_func{:});
    end
    toc
    tic
    for test_id=1:100
        Phi_func(x_initial_cell{:},param_initial,other_input_Phi_func{:});
    end
    toc

    temp="XXX"
    