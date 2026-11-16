function [obj_val,grad]=obj_nfxp_func(...
    param,Phi_func,obj_func,...
other_vars_Phi_func,other_vars_obj_func,spec,varargin)

global x_initial_cell_global count_Phi_eval
global x_sol_cell
%param'

x_initial0_cell=varargin;

if isempty(x_initial_cell_global)==0
    x_initial_cell=x_initial_cell_global;%hot start
else
    x_initial_cell=x_initial0_cell;
end

hot_start_spec=1;
if hot_start_spec==0 %Always start from x_initial0_cell, which is the initial guess specified by the user.
    x_initial_cell=x_initial0_cell;%%%%%
end

n_vars=size(x_initial_cell,2);

if isempty(spec)==1
    spec=[];
    spec.TOL=1e-12;
    spec.ITER_MAX=2000;
    spec.Anderson_acceleration=1;
end

if isfield(spec,'nfxp_Newton_spec')==0
    spec.TOL=1e-12;%%%%%
    spec.ITER_MAX=2000;%%%%%
    spec.Anderson_acceleration=1;%%%%

    [output_spectral,other_vars,iter_info]=...
            spectral_func(Phi_func,spec,...
            x_initial_cell,param,other_vars_Phi_func{:});

    if iter_info.FLAG_ERROR==1 | iter_info.feval>=iter_info.ITER_MAX
        %temp="Not converge"
        %iter_info.feval
        %iter_info.DIST_table(end,:)
    end

    x_sol_cell=output_spectral;
    %iter_info.feval
    count_Phi_eval=count_Phi_eval+iter_info.feval;

else % Use Newton iteration to solve for fixed point
    spec_gmres=[];

    if isfield(spec,'Newton_update_func')==1
        spec_gmres.Newton_update_func=spec.Newton_update_func;
    else
        spec_gmres.Newton_update_func=[];
    end

    n_vars=size(x_initial_cell(:),1);
    fixed_point_iter_spec=1;
    [x_sol_cell,iter_info_x,other_vars,flag_vec,relres_vec,iter_cell,resvec_cell]=...
            Newton_func(x_initial0_cell{:},Phi_func,fixed_point_iter_spec,...
                    param,other_vars_Phi_func{:},spec_gmres,@L2_norm_func,n_vars);
end

    x_initial_cell_global=x_sol_cell;


%% Compute objective value (and gradient if needed)
    krylov_spec=0;
    if nargout==1
        obj_val=obj_func(param,x_sol_cell{:},other_vars_obj_func{:});

        %obj_val
        if isnan(obj_val)==1|isinf(obj_val)==1 | iter_info.n_iter==spec.ITER_MAX
            obj_val=1e+10;
            x_initial_cell_global=[];
        end

    elseif nargout>=2 & isfield(spec,'compute_z_y_func')==1
        %% Analytically compute gradient
        [obj_val,dobj_dparam_direct,dobj_dX_cell]=obj_func(param,x_sol_cell{:},other_vars_obj_func{:});

        %obj_val
        if isnan(obj_val)==1|isinf(obj_val)==1 | iter_info.n_iter==spec.ITER_MAX
            obj_val=1e+10;
            x_initial_cell_global=[];
        end

        compute_z_y_func=spec.compute_z_y_func;
        [y_cell]=compute_z_y_func(x_sol_cell{:},param,other_vars_Phi_func{:});
        
        n_param=size(param(:),1);
        grad=reshape(dobj_dparam_direct,1,n_param);
        for id=1:n_vars
          grad=grad+...
              reshape(dobj_dX_cell{id},[],1)'*...
              reshape(y_cell{id},[],n_param);
        end


    elseif nargout>=2 & krylov_spec==1 % Compute gradient by numerical derivatives
      Richardson_spec=0;
      fixed_point_iter_spec=1;%%%
      n_vars=size(x_initial_cell,1);
      n_param=size(param,1);

      dPhi_dparam_k_cell=compute_dPhi_dparam_direct_func(...
        x_sol_cell,param,Phi_func,other_vars_Phi_func,Richardson_spec);

      for i=1:n_vars
        if fixed_point_iter_spec==1
            dG_dparam_k_cell{i}=-dPhi_dparam_k_cell{i};
        else% fixed_point_iter_spec==0
            dG_dparam_k_cell{i}=dPhi_dparam_k_cell{i};
        end

        dG_dparam_k_cell{i}=reshape(dG_dparam_k_cell{i},[],n_param);
        dG_dparam_k_cell_negative{i}=-dG_dparam_k_cell{i};

        y_cell_init{i}=dG_dparam_k_cell_negative{i}*0;%%%
      end

      spec_gmres=[];
      spec_gmres.ITER_MAX=200;
      [y_cell,flag_vec_y,relres_vec_y,iter_cell_y,resvec_cell_y]=...
            solve_linear_eq_JF_cell_func(dG_dparam_k_cell_negative,y_cell_init,spec_gmres,Phi_func,n_vars,...
            x_sol_cell{:},param,other_vars_Phi_func{:},...
                fixed_point_iter_spec);

      for i=1:n_vars
        size_x{i}=size(x_sol_cell{i});
      end
  
      dobj_dparam_indirect=compute_dobj_dparam_indirect_func(...
          param,x_sol_cell,y_cell,size_x,obj_func,other_vars_obj_func,Richardson_spec,n_vars);

      dobj_dparam_direct=compute_dobj_dparam_direct_func(...
          param,x_sol_cell,obj_func,other_vars_obj_func,Richardson_spec);

      grad=dobj_dparam_direct+dobj_dparam_indirect;

    else

        grad=[];
    end% Compute gradient by solving linear equations

    
end
