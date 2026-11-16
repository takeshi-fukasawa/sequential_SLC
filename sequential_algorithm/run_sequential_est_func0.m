function [param_SLC,iter_info_SLC,X_sol_SLC,...
    iter_info_x,other_vars,speed_test_info,results_nfxp,results_EPL]=...
    run_sequential_est_func0(param_initial,varargin)

    
    warning('off', 'all');

    global count_obj_eval
    global G_scaling_cell x_sol_cell
    global positive_param_nfxp_spec
    global x_initial_cell_global
    
    n_vars_x=varargin{end};
    n_vars=n_vars_x;
    
    for m=1:n_vars
        x_initial0_cell{m}=varargin{m};
        G_scaling_cell{m}=1;
    end

    c=0;
    c=c+1;spec_gmres=varargin{end-c};
    c=c+1;spec_sequential=varargin{end-c};
    c=c+1;opt_options=varargin{end-c};
    c=c+1;rho_val0=varargin{end-c};
    c=c+1;other_input_obj_func=varargin{end-c};
    c=c+1;other_input_Phi_func=varargin{end-c};
    c=c+1;middle_vars_spec_cell=varargin{end-c};
    c=c+1;obj_func=varargin{end-c};
    c=c+1;Phi_func=varargin{end-c};

    central_difference_spec=1;
    fixed_point_iter_spec=1; %% Currently, Phi_func should be a fixed point mapping. 


if  isfield(spec_sequential,'Richardson_spec')==0
    Richardson_spec=0;% If 1, use Richardson extrapolation. If 0, do not use it.
else
    Richardson_spec=spec_sequential.Richardson_spec;
end

    if isfield(spec_sequential,"SLC_positive_param_update_spec")==0
        spec_sequential.SLC_positive_param_update_spec=0;
    end


if  isfield(spec_sequential,'line_search_spec')==0
    spec_sequential.line_search_spec=0;
end

if  isfield(spec_sequential,'initial_x_solve_spec')==0
    spec_sequential.initial_x_solve_spec=1;
end

%if  isfield(spec_sequential,'no_spectral_spec')==1
%    if spec_sequential.no_spectral_spec==1
%        spec_sequential.spectral_coef_spec=0;
%    end
%end

if isfield(spec_sequential,'spectral_spec')==0
    spectral_spec=0;
else
    spectral_spec=spec_sequential.spectral_spec;
end

    if isfield(spec_sequential,'x_min_cell')==1
       param_min=spec_sequential.x_min_cell{1};
       param_initial=max(param_initial,param_min);
    else
         param_min=[];
    end

    if isfield(spec_sequential,'x_max_cell')==1
       param_max=spec_sequential.x_max_cell{1};
       param_initial=min(param_initial,param_max);
    else
         param_max=[];
    end

if isfield(spec_sequential,'Newton_update_func')==1
    spec_gmres.Newton_update_func=spec_sequential.Newton_update_func;
end

if isfield(spec_sequential,'run_EPL_spec')==0
    spec_sequential.run_EPL_spec=0;
end

%% Solve x_initial
spec=[];
spec.TOL=1e-12;

    spec.ITER_MAX=2000;%%%%%
    spec.Anderson_acceleration=1;
    %spec.update_spec=0;
    spec.DEBUG=1;
    
        if isempty(other_input_Phi_func)==1
            [x_sol_cell,~,iter_info_x]=...
            spectral_func(Phi_func,spec,{x_initial0_cell{:}},param_initial);
        else
            [x_sol_cell,~,iter_info_x]=...
            spectral_func(Phi_func,spec,{x_initial0_cell{:}},param_initial,other_input_Phi_func{:});
        end

        if 1==0 %% Try Newton iteration to solve for fixed point
            if isempty(other_input_Phi_func)==1
                [x_sol_cell,iter_info_x,other_vars,flag_vec,relres_vec,iter_cell,resvec_cell]=...
                    Newton_func(x_initial0_cell{:},Phi_func,fixed_point_iter_spec,...
                    param_initial,spec_gmres,@L2_norm_func,n_vars);
            else
                [x_sol_cell,iter_info_x_Newton,other_vars,flag_vec,relres_vec,iter_cell,resvec_cell]=...
                    Newton_func(x_initial0_cell{:},Phi_func,fixed_point_iter_spec,...
                    param_initial,other_input_Phi_func{:},spec_gmres,@L2_norm_func,n_vars);
            end
        end
    x_initial_cell=x_sol_cell;
    
    %% Use original initial value

    if spec_sequential.initial_x_solve_spec==0|iter_info_x.FLAG_ERROR==1
        x_initial_cell=x_initial0_cell;
    end

    %% Measure comp time for evaluating obj_func
     other_input_obj_func_test=other_input_obj_func;

    if isempty(middle_vars_spec_cell)==0
        middle_vars_func=middle_vars_spec_cell.middle_vars_func;
        middle_vars_cell=middle_vars_func(x_initial_cell{:},middle_vars_spec_cell.other_input_middle_vars_func{:});
        other_input_obj_func_test{1}=middle_vars_cell;
    end


    n_iter_test=100;%100;
    tStart_test = tic; 
    for test_id=1:n_iter_test
        [obj_test,~]=obj_func(param_initial,x_initial_cell{:},other_input_obj_func_test{:});
    end
    t_test=toc(tStart_test)
    speed_test_info.comp_time_eval_obj_func=t_test/n_iter_test;
    
    speed_test_info.comp_time_eval_Phi_func=iter_info_x.t_cpu/iter_info_x.n_iter;
    

    if isfield(spec_sequential,'analytical_spec')==0
        spec_sequential.analytical_spec=0;
    end
    
    %%%%%%%%%%%%%%%%%
    %% NFXP
    "NFXP"
    if isfield(spec_sequential,"spec_inner_iter")==0
        spec_inner_iter=[];
        spec_inner_iter.Anderson_acceleration=1;

        if isfield(spec_sequential,"compute_z_y_func")==1 & spec_sequential.analytical_spec==1
            spec_inner_iter.compute_z_y_func=spec_sequential.compute_z_y_func;
        end
    end

    results_nfxp=[];
    if isfield(spec_sequential,"run_nfxp")==0
        spec_sequential.run_nfxp=0;
    end
    if spec_sequential.run_nfxp==1

        opt_options_nfxp=optimset(...
            'FinDiffType','central','MaxIter',200,...
                    'Display','iter','TolFun',1e-6,...
                        'TolX',1e-6,'FunValCheck','off',MaxFunEvals = 2000);

        if isfield(spec_sequential,"nfxp_bhhh_spec")==1
            if spec_sequential.nfxp_bhhh_spec==1
                spec_nfxp.bhhh_spec=1;
                opt_options_nfxp.TolFun=1e-1;
                opt_options_nfxp.TolX=1e-1;
            else
                spec_nfxp.bhhh_spec=0;
            end
        else
            spec_nfxp.bhhh_spec=0;
        end

        if spec_sequential.analytical_spec==1
            opt_options_nfxp.GradObj = 'on';
        else
            opt_options_nfxp.GradObj = 'off';
        end

        if isfield(spec_sequential,"constraint_func")==0
            spec_nfxp.constraint_func=[];
        else
            spec_nfxp.constraint_func=spec_sequential.constraint_func;
        end
        if isfield(spec_sequential,"x_min_cell")==1
            spec_nfxp.param_min=spec_sequential.x_min_cell{1};
        end
        if isfield(spec_sequential,"x_max_cell")==1
            spec_nfxp.param_max=spec_sequential.x_max_cell{1};
        end

        if isfield(spec_inner_iter,"compute_z_y_func")==1
            spec_inner_iter.compute_z_y_func=spec_sequential.compute_z_y_func;
        end
        if isfield(spec_inner_iter,"Newton_update_func")==1
            spec_inner_iter.Newton_update_func=spec_sequential.Newton_update_func;
        end
            
        if isfield(spec_sequential,"nfxp_Newton_spec")==1
            spec_inner_iter.nfxp_Newton_spec=spec_sequential.nfxp_Newton_spec;
        end

        if spec_sequential.SLC_positive_param_update_spec==1
            positive_param_nfxp_spec=1;
        else
            positive_param_nfxp_spec=0;
        end

        [param_nfxp,iter_info_nfxp]=...
            run_nfxp_func(param_initial,spec_inner_iter,spec_nfxp,x_initial0_cell{:},...
            opt_options_nfxp,obj_func,Phi_func,...
            other_input_obj_func,other_input_Phi_func,n_vars);
        
        results_nfxp.param_nfxp=param_nfxp;
        results_nfxp.iter_info_nfxp=iter_info_nfxp;
    else
        results_nfxp=[];
    end


    %%%%% SLC %%%%%
        "SLC"
        global count_G_eval
        global y_cell_init z_cell_init
        
        positive_param_nfxp_spec=0;
        
        count_obj_eval=0;
        count_G_eval=0;
        
        opt_options_SLC = optimset('FinDiffType','central','MaxIter',200,...
            'Display','off','TolFun',1e-7,...
            'TolX',1e-7,'FunValCheck','off',MaxFunEvals = 2000);

        if spec_sequential.analytical_spec==1
            opt_options_SLC.GradObj = 'on';
        else          
            opt_options_SLC.GradObj = 'off';
        end

        
            spec_SLC_algorithm=spec_sequential;
            spec_SLC_algorithm.opt_options=opt_options_SLC;
            spec_SLC_algorithm.spec_gmres.TOL=1e-8;
            spec_SLC_algorithm.spec_gmres.ITER_MAX=500;
            spec_SLC_algorithm.common_spectral_coef_spec=1;
            spec_SLC_algorithm.obj_func=obj_func;
            spec_SLC_algorithm.Phi_func=Phi_func;
            spec_SLC_algorithm.other_input_obj_func=other_input_obj_func;
            spec_SLC_algorithm.other_input_Phi_func=other_input_Phi_func;
            spec_SLC_algorithm.central_difference_spec=central_difference_spec;
            spec_SLC_algorithm.fixed_point_iter_spec=fixed_point_iter_spec;
    
            spec_SLC_algorithm.Richardson_spec=Richardson_spec;
            spec_SLC_algorithm.param_min=param_min;
            spec_SLC_algorithm.param_max=param_max;
        
            spec_SLC_algorithm.x_min_cell=cell(1,size(x_initial_cell,2)+1);
            spec_SLC_algorithm.x_max_cell=cell(1,size(x_initial_cell,2)+1);
            if isempty(param_min)==0
                spec_SLC_algorithm.x_min_cell{1}=param_min;
            end
            if isempty(param_max)==0
                spec_SLC_algorithm.x_max_cell{1}=param_max;
            end

            if isfield(spec_sequential,"compute_z_y_func")==1
                if spec_sequential.analytical_spec==1
                    spec_SLC_algorithm.compute_z_y_func=spec_sequential.compute_z_y_func;
                else
                    spec_SLC_algorithm.compute_z_y_func=[];
                end
            else
                spec_SLC_algorithm.compute_z_y_func=[];
            end



            %if iter_info_x.FLAG_ERROR==1 % If x_initial is not solved accurately, start from original x_initial
            %    x_initial_cell=x_initial0_cell;
            %end
            
            if 1==0 & spec_sequential.run_nfxp==1 %% test: Start from NFXP solution
                x_initial_cell=x_initial_cell_global;%%%%
                param_initial=param_nfxp*1.0;%%%%%
            end
            
            if isfield(spec_sequential,"TOL")==0
                spec_SLC_algorithm.TOL=1e-6;%%%%
            end

            spec_SLC_algorithm.ITER_MAX=50;%100;%30;
            spec_SLC_algorithm.line_search_spec=spec_sequential.line_search_spec;
            %spec_SLC_algorithm.norm_spec=0;%%%%%

            if spectral_spec==0
                spec_SLC_algorithm.update_spec=0;%Use original SLC algorithm (Not introduce spectral). 
            elseif spectral_spec==1
                spec_SLC_algorithm.update_spec=[];%Use spectral SLC algorithm.
            elseif spectral_spec==2
                spec_SLC_algorithm.Anderson_acceleration=1;%%%%%% Use Anderson acceleration for SLC. 
            end
        


         y_cell_init=[];z_cell_init=[];

         %x_initial_cell{1}=x_initial_cell{1}*0;
         %param_initial=param_initial*1;
         X_SLC_initial_cell={param_initial,x_initial_cell{:}};

         %spectral_func_for_constrained_opt
         [X_sol_cell,other_vars,iter_info_SLC]=...
                spectral_func_for_constrained_opt(@SLC_update_func,spec_SLC_algorithm,X_SLC_initial_cell,...
                Phi_func,obj_func,...
                other_input_Phi_func,other_input_obj_func,spec_SLC_algorithm,n_vars);
        param_SLC=X_sol_cell{1};

        X_sol_SLC=X_sol_cell;

        iter_info_SLC.count_obj_eval=count_obj_eval;
        iter_info_SLC.count_G_eval=count_G_eval;

        iter_info_SLC.param_initial=X_SLC_initial_cell{1};

        if spec_sequential.initial_x_solve_spec==1
            iter_info_SLC.t_cpu_sequential=iter_info_SLC.t_cpu;
            iter_info_SLC.t_cpu=iter_info_SLC.t_cpu+iter_info_x.t_cpu;
            iter_info_SLC.t_cpu_sol_x=iter_info_x.t_cpu;

            iter_info_SLC.count_G_eval=iter_info_SLC.count_G_eval+iter_info_x.feval;            
        end

        [obj_val_SLC,~]=obj_func(param_SLC,X_sol_cell{2:end},other_input_obj_func{:});

        iter_info_SLC.obj_val=obj_val_SLC;

    if spec_sequential.run_EPL_spec==1 % Run EPL algorithm
        "EPL"

        spec_EPL_algorithm=spec_SLC_algorithm;
        spec_EPL_algorithm.ITER_MAX=50;
        spec_EPL_algorithm.TOL=1e-6;%%%%

        opt_options_EPL = optimset('FinDiffType','central','MaxIter',200,...
            'Display','off','TolFun',1e-6,...
            'TolX',1e-6,'FunValCheck','off',MaxFunEvals = 2000);
        spec_EPL_algorithm.opt_options=opt_options_EPL;

        y_cell_init=[];z_cell_init=[];

         %x_initial_cell{1}=x_initial_cell{1}*0;
         %param_initial=param_initial*1;
         X_EPL_initial_cell={param_initial,x_initial_cell{:}};

         [X_sol_cell,other_vars,iter_info_EPL]=...
                spectral_func_for_constrained_opt(@EPL_update_func,spec_EPL_algorithm,X_EPL_initial_cell,...
                Phi_func,obj_func,...
                other_input_Phi_func,other_input_obj_func,spec_EPL_algorithm,n_vars);
        param_EPL=X_sol_cell{1};

        X_sol_EPL=X_sol_cell;

        iter_info_EPL.count_obj_eval=count_obj_eval;
        iter_info_EPL.count_G_eval=count_G_eval;

        iter_info_EPL.param_initial=X_EPL_initial_cell{1};

        [obj_val_EPL,~]=obj_func(param_EPL,X_sol_cell{2:end},other_input_obj_func{:});

        iter_info_EPL.obj_val=obj_val_EPL;

        results_EPL.iter_info_EPL=iter_info_EPL;
        results_EPL.param_EPL=param_EPL;
        
    else
        results_EPL=[];
    end

end
