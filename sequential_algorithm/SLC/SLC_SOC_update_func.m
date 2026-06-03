function [out,other_vars]=...
    SLC_SOC_update_func(param_k,varargin)


    global k_iter output_opt count_obj_eval
    global param_initial_temp
    global count_G_eval
    global y_cell_init z_cell_init
    global G_scaling_cell
    
    rho_val=0.0;% Weight used in the EPL objective function; Currently set to 0

    
    param_k'

    %% Step 0. Preliminary
    if isempty(k_iter)==1
        k_iter=1;
     end

    n_vars=varargin{end};
    for m=1:n_vars
        x_k_cell{m}=varargin{m};
        cell_size{m}=size(x_k_cell{m});
    end

    c=0;
    c=c+1;spec_epl_algorithm=varargin{end-c};
    c=c+1;other_input_obj_func0=varargin{end-c};
    c=c+1;other_input_Phi_func0=varargin{end-c};
    c=c+1;obj_func_given_x=varargin{end-c};
    c=c+1;Phi_func=varargin{end-c};

    spec_gmres=spec_epl_algorithm.spec_gmres;
    Richardson_spec=spec_epl_algorithm.Richardson_spec;
    central_difference_spec=spec_epl_algorithm.central_difference_spec;
    
    param_min=spec_epl_algorithm.param_min;
    param_max=spec_epl_algorithm.param_max;
    opt_options=spec_epl_algorithm.opt_options;
    fixed_point_iter_spec=spec_epl_algorithm.fixed_point_iter_spec;

    if isfield(spec_epl_algorithm,'constraint_func')==1
        constraint_func=spec_epl_algorithm.constraint_func;
    else
        constraint_func=[];
    end

    if isempty(param_min)==0
        param_k=max(param_k,param_min);
    end
    if isempty(param_max)==0
        param_k=min(param_k,param_max);
    end
    


    param_initial_temp=param_k;

    n_param=size(param_k(:),1);

    if isempty(z_cell_init)==0
        if isnan(sum(z_cell_init{1}(:)))==1
            z_cell_init=[];y_cell_init=[];
        end
    end
    
    %% Step 1: Update x



        dPhi_dparam_k_cell=compute_dPhi_dparam_direct_func(...
            x_k_cell,param_k,Phi_func,other_input_Phi_func0,Richardson_spec);

        for i=1:n_vars
            if fixed_point_iter_spec==1
                dG_dparam_k_cell{i}=-dPhi_dparam_k_cell{i};
            else% fixed_point_iter_spec==0
                dG_dparam_k_cell{i}=dPhi_dparam_k_cell{i};
            end

        dG_dparam_k_cell{i}=reshape(dG_dparam_k_cell{i},[],n_param);

        dG_dparam_k_cell{i}=dG_dparam_k_cell{i}./G_scaling_cell{i};
        dG_dparam_k_cell_negative{i}=-dG_dparam_k_cell{i};
        end

        dG_dparam_k_mat=[];%%[]*n_param
        for m=1:n_vars
            dG_dparam_k_mat=[dG_dparam_k_mat;...
                dG_dparam_k_cell{m}];
        end


        spec_gmres.ITER_MAX=300;%%%
        [y_k_cell,flag_vec_y,relres_vec_y,iter_cell_y,resvec_cell_y]=...
                solve_linear_eq_JF_cell_func(dG_dparam_k_cell_negative,y_cell_init,spec_gmres,Phi_func,n_vars,...
               x_k_cell{:},param_k,other_input_Phi_func0{:},...
                    fixed_point_iter_spec);
    
    
        for i=1:size(iter_cell_y,2)
            count_G_eval=count_G_eval+prod(iter_cell_y{i})*2;
        end
        %y_k_cell_temp=y_k_cell;
        %z_k_cell_temp=z_k_cell;
        %G_k_cell_temp=G_k_cell;   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SOC_spec=0;%% 0: without SOC; 1,2,3,...: with SOC (Second-order correction); 

    if k_iter>=10
        SOC_spec=1;
    else
        SOC_spec=0;
    end

    SOC_spec=0;%%%

    for iter_SOC=1:SOC_spec+1
        %% Compute G_k_cell

        G_k_cell=G_func(Phi_func,x_k_cell,fixed_point_iter_spec,param_k,other_input_Phi_func0{:});

        for m=1:n_vars
            G_k_cell_negative{m}=-G_k_cell{m};
        end

        count_G_eval=count_G_eval+1;

        spec_gmres.ITER_MAX=300;
        [z_k_cell,flag_vec_z,relres_vec_z,iter_cell_z,resvec_cell_z]=...
           solve_linear_eq_JF_cell_func(G_k_cell_negative,z_cell_init,spec_gmres,Phi_func,n_vars,...
           x_k_cell{:},param_k,other_input_Phi_func0{:},...
                fixed_point_iter_spec);

        for i=1:size(iter_cell_z,2)
            count_G_eval=count_G_eval+prod(iter_cell_z{i});
        end

    
        z_cell_init=z_k_cell;%% Used in the next iteration
        y_cell_init=y_k_cell;%% Used in the next iteration

    
        %% Step 3: Update param

        x_k_temp_cell=x_k_cell;
        y_k_temp_cell=y_k_cell;
        z_k_temp_cell=z_k_cell;


        %opt_options = optimset('GradObj','off',...
            %'FinDiffType','central','FinDiffRelStep',1e-8,...
            %'TolFun',1e-8);


            if isempty(param_min)==1 & isempty(param_max)==1
                try
                    [param_updated,fval,exitflag,output_opt]=...
                    fminunc(@obj_func_SLC,param_k*1.0,opt_options,...
                    param_k,obj_func_given_x,...
                    x_k_temp_cell,y_k_temp_cell,z_k_temp_cell,other_input_obj_func0,...
                    Phi_func,other_input_Phi_func0,fixed_point_iter_spec,G_k_cell,rho_val,constraint_func);
                catch
                    param_updated=param_k*NaN;
                    output_opt.funcCount=NaN;
                end
            else%with upper/lower bounds
                try
                    [param_updated,fval,exitflag,output_opt]=...
                        fmincon(@obj_func_SLC,param_k*1.0,[],[],[],[],...
                        param_min,param_max,[],...
                        opt_options,...
                        param_k,obj_func_given_x,...
                        x_k_temp_cell,y_k_temp_cell,z_k_temp_cell,other_input_obj_func0,...
                        Phi_func,other_input_Phi_func0,fixed_point_iter_spec,G_k_cell,rho_val,constraint_func);
                catch
                    param_updated=param_k*NaN;
                    output_opt.funcCount=NaN;
                end
            end

           % if max(abs(param_updated-param_k))<1e-5
           %     param_updated=param_k;
           % end
           
           if isfield(opt_options,'GradObj')==1 && strcmp(opt_options.GradObj,'on')==1
                count_obj_eval=count_obj_eval+(output_opt.funcCount);
           else % Numerical derivative based 
                count_obj_eval=count_obj_eval+(output_opt.funcCount)*(2*n_param+1);
           end


        %output_opt
        %if exitflag~=1 && exitflag~=2
        %    output_opt
        %end


        %param_updated=param_k;%%%%

        %% Step 2: Update x
        for m=1:n_vars
            x_updated_cell{m}=x_k_cell{m}+z_k_cell{m}+...
                reshape(y_k_cell{m}*(param_updated-param_k),size(x_k_cell{m}));
        end

        if iter_SOC==1 & SOC_spec>=1
            param_k=param_updated;
            x_k_cell=x_updated_cell;
        end

    end% for loop wrt iter_SOC


   out={param_updated,x_updated_cell{:}};

    k_iter=k_iter+1;

    other_vars=[];

end
