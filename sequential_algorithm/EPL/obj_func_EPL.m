function [obj_val,grad]=obj_func_EPL(param,...
    param_initial,obj_fun,x_cell,Jac,other_input_obj_func,...
    Phi_func,other_input_Phi_func0,fixed_point_iter_spec,spec_gmres,compute_dG_dparam_func)

    global z_cell_init count_G_eval

        n_vars=size(x_cell,2);
        n_param=size(param(:),1);

      for i=1:n_vars
        cell_size{i}=size(x_cell{i});
      end

     G_cell=G_func(Phi_func,x_cell,fixed_point_iter_spec,param,other_input_Phi_func0{:});

     count_G_eval=count_G_eval+1;
     
        for m=1:n_vars
            G_cell_negative{m}=-G_cell{m};
        end

    if nargout==1
          [z_cell,flag_vec_z,relres_vec_z,iter_cell_z,resvec_cell_z]=...
              solve_linear_eq_JF_cell_func(G_cell_negative,z_cell_init,spec_gmres,Phi_func,n_vars,...
              x_cell{:},param,other_input_Phi_func0{:},...
                      fixed_point_iter_spec);

          for i=1:size(iter_cell_z,2)
              count_G_eval=count_G_eval+prod(iter_cell_z{i});
          end

          z_cell_init=z_cell;%% Used in the next iteration
    else% nargout==2
          G_vec=cell_to_vec_func(G_cell);
       
        z_vec=-Jac\G_vec;

        z_cell=vec_to_cell_func1(z_vec,cell_size);
    end

    for m=1:n_vars
      x_cell_for_obj_func{m}=x_cell{m}+z_cell{m};
    end%for m=1:n_vars

    if nargout==1
          if isempty(other_input_obj_func)==1
            [obj_val]=obj_fun(param,x_cell_for_obj_func{:});
          else
            [obj_val]=obj_fun(param,x_cell_for_obj_func{:},other_input_obj_func{:});
          end

    else%nargout==2 

        if isempty(other_input_obj_func)==1
          [obj_val,dobj_dparam_direct,dobj_dX_cell]=obj_fun(param,x_cell_for_obj_func{:});
        else
          [obj_val,dobj_dparam_direct,dobj_dX_cell]=obj_fun(param,x_cell_for_obj_func{:},other_input_obj_func{:});
        end

        dG_dparam=compute_dG_dparam_func(x_cell_for_obj_func{:},param,other_input_Phi_func0{:});   
        dX_dparam=-Jac\dG_dparam;

        dobj_dX_vec=[];
        for id=1:n_vars
           dobj_dX_vec=[dobj_dX_vec;dobj_dX_cell{id}(:)];
       end

       grad=reshape(dobj_dparam_direct,1,n_param)+...
            reshape(dobj_dX_vec,1,[])*dX_dparam;%1*n_param

    end%nargout==1 or not

    %obj_val
end
