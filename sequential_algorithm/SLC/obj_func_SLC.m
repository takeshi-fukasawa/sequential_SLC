function [obj_val,grad]=obj_func_SLC(param,...
    param_initial,obj_fun,x_cell,y_cell,z_cell,other_input_obj_func,Phi_func,other_input_Phi_func0,fixed_point_iter_spec,G_k_cell_0,rho_val,constraint_func)

    global deviation G_k_cell_0_global

    G_k_cell_0_global=G_k_cell_0;

        n_vars=size(x_cell,2);
        n_param=size(param(:),1);

        for m=1:n_vars
            x_cell_for_obj_func{m}=x_cell{m}+z_cell{m}+...
                reshape(reshape(y_cell{m},[],n_param)*(param-param_initial),size(x_cell{m}));
        end%for m=1:n_vars

      if nargout==1
          if isempty(other_input_obj_func)==1
            [obj_val]=obj_fun(param,x_cell_for_obj_func{:});
          else
            [obj_val]=obj_fun(param,x_cell_for_obj_func{:},other_input_obj_func{:});
          end
      else% Compute gradient
        if isempty(other_input_obj_func)==1
          [obj_val,dobj_dparam_direct,dobj_dX_cell]=obj_fun(param,x_cell_for_obj_func{:});
        else
          [obj_val,dobj_dparam_direct,dobj_dX_cell]=obj_fun(param,x_cell_for_obj_func{:},other_input_obj_func{:});
        end
        grad=reshape(dobj_dparam_direct,1,n_param);
        for id=1:n_vars
          grad=grad+...
              reshape(dobj_dX_cell{id},[],1)'*...
              reshape(y_cell{id},[],n_param);
        end

      end

      %rho_val=100;%0.1;%%%
      %rho_val=10;%0.1;%%%
      %rho_val=0;%0.1;%%%
      

      rho_val_proximal=0.0;
      obj_val=obj_val+rho_val_proximal*sum((param(:)-param_initial(:)).^2); % Now work well...

      if rho_val>0

        if isempty(constraint_func)==0
          G_k_cell=constraint_func(x_cell_for_obj_func{:},param,other_input_Phi_func0{:});
        else
          G_k_cell=G_func(Phi_func,x_cell_for_obj_func,fixed_point_iter_spec,param,other_input_Phi_func0{:});
        end

        %rho_val=0;%%%
        for i=1:size(G_k_cell,2)
          obj_val=obj_val+rho_val*sum(G_k_cell{i}(:).^2);
          
          %deviation=G_k_cell_0{i}(:)-dG_dparam_cell_0{i}*(param-param_initial);
          %deviation=z_cell{i}(:)+y_cell{i}*(param-param_initial);
          
          %obj_val=obj_val+rho_val*sum(deviation(:).^2); % Now work well...
        end

        %obj_val=obj_val+rho_val*sum(G_k_cell{1}(:).^2);
        %obj_val=obj_val+rho_val*sum(G_k_cell{1}(1).^2);
        

      end

end
