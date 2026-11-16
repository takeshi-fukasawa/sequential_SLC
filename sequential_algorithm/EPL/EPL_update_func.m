function [out,other_vars]=...
    EPL_update_func(param_k,varargin)


    global k_iter output_opt count_obj_eval
    global param_initial_temp
    global count_G_eval
    global y_cell_init z_cell_init
    global G_scaling_cell
    
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
    c=c+1;spec_EPL_algorithm=varargin{end-c};
    c=c+1;other_input_obj_func0=varargin{end-c};
    c=c+1;other_input_Phi_func0=varargin{end-c};
    c=c+1;obj_func_given_x=varargin{end-c};
    c=c+1;Phi_func=varargin{end-c};

    spec_gmres=spec_EPL_algorithm.spec_gmres;

    param_min=spec_EPL_algorithm.param_min;63333
    param_max=spec_EPL_algorithm.param_max;
    opt_options=spec_EPL_algorithm.opt_options;
    fixed_point_iter_spec=spec_EPL_algorithm.fixed_point_iter_spec;

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

    if isfield(spec_EPL_algorithm,'compute_Jac_func')==0
        spec_EPL_algorithm.compute_Jac_func=[];
    end
    
    if isempty(spec_EPL_algorithm.compute_Jac_func)==1 
        dJac=[];
        opt_options.GradObj='off';
        compute_dG_dparam_func=[];
    else
        compute_Jac_func=spec_EPL_algorithm.compute_Jac_func;
        compute_dG_dparam_func=spec_EPL_algorithm.compute_dG_dparam_func;
        Jac=compute_Jac_func(x_k_cell{:},param_k,other_input_Phi_func0{:});
        
         opt_options.GradObj='on';
         dJac=decomposition(Jac);

         compute_dG_dparam_func=spec_EPL_algorithm.compute_dG_dparam_func;
    end
        
   %% Step 1: Update param

   param_k'

       x_k_temp_cell=x_k_cell;


    %opt_options = optimset('GradObj','off',...
        %'FinDiffType','central','FinDiffRelStep',1e-8,...
        %'TolFun',1e-8);

        %opt_options.GradObj='off';

            if isempty(param_min)==1 & isempty(param_max)==1
                try
                    [param_updated,fval,exitflag,output_opt]=...
                    fminunc(@obj_func_EPL,param_k*1.0,opt_options,...
                    param_k,obj_func_given_x,...
                    x_k_cell,dJac,other_input_obj_func0,...
                    Phi_func,other_input_Phi_func0,fixed_point_iter_spec,spec_gmres,compute_dG_dparam_func);
                catch
                    param_updated=param_k*NaN;
                    output_opt.funcCount=NaN;
                end
            else%with upper/lower bounds
                try
                    [param_updated,fval,exitflag,output_opt]=...
                        fmincon(@obj_func_EPL,param_k*1.0,[],[],[],[],...
                        param_min,param_max,[],...
                        opt_options,...
                        param_k,obj_func_given_x,...
                        x_k_cell,dJac,other_input_obj_func0,...
                        Phi_func,other_input_Phi_func0,fixed_point_iter_spec,spec_gmres,compute_dG_dparam_func);
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
                count_obj_eval=count_obj_eval+(output_opt.funcCount);
           end

        %output_opt
        %if exitflag~=1 && exitflag~=2
        %    output_opt
        %end

        if spec_EPL_algorithm.SLC_positive_param_update_spec==1
            param_updated=abs(param_updated);
        end

    %param_updated=param_k;%%%%

    %% Step 2: Update x

    G_updated_cell=G_func(Phi_func,x_k_cell,fixed_point_iter_spec,param_updated,other_input_Phi_func0{:});

    for m=1:n_vars
        G_updated_cell_negative{m}=-G_updated_cell{m};
    end

    count_G_eval=count_G_eval+1;

    if isempty(spec_EPL_algorithm.compute_Jac_func)==1 

        [z_updated_cell,flag_vec_z,relres_vec_z,iter_cell_z,resvec_cell_z]=...
            solve_linear_eq_JF_cell_func(G_updated_cell_negative,z_cell_init,spec_gmres,Phi_func,n_vars,...
            x_k_cell{:},param_k,other_input_Phi_func0{:},...
                    fixed_point_iter_spec);

        for i=1:size(iter_cell_z,2)
            count_G_eval=count_G_eval+prod(iter_cell_z{i});
        end

        z_cell_init=z_updated_cell;%% Used in the next iteration

    elseif isempty(spec_EPL_algorithm.compute_z_y_func)==0
        G_updated_vec=[];
        for m=1:n_vars
            G_updated_vec=[G_updated_vec;G_updated_cell{m}(:)];
       end

         z_updated_vec=-dJac\G_updated_vec;
         z_updated_cell=vec_to_cell_func1(z_updated_vec,cell_size);
            
    end

    for m=1:n_vars
        x_updated_cell{m}=x_k_cell{m}+z_updated_cell{m};
    end


   out={param_updated,x_updated_cell{:}};

    k_iter=k_iter+1;

    other_vars=[];

end
