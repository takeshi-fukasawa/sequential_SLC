%%%%%%%%%%%%%%%
function dobj_dparam_indirect=compute_dobj_dparam_indirect_func(param,x_cell,dx_dparam_cell,size_x,obj_func_given_x,other_input_obj_func0,Richardson_spec,n_vars)

    dobj_dparam_indirect=zeros(size(param));

    n_param=size(param(:),1);
    eps_val=eps^(1/3);

    for m=1:n_vars
        dx_dparam_reshaped_cell{m}=reshape(dx_dparam_cell{m},size(x_cell{m}(:),1),size(param(:),1));
    end
    for i=1:n_param
        norm_dx_dparam_i=zeros(n_vars,1);
        for m=1:n_vars
            norm_dx_dparam_i(m)=max(abs(dx_dparam_reshaped_cell{m}(:,i)));
            dx_dparam_i_cell{m}=reshape(dx_dparam_reshaped_cell{m}(:,i),size_x{m});
        end

        eps_val_scaled_i=eps_val./max(max(norm_dx_dparam_i),1e-8);%scaling
        %eps_val_scaled_i=eps_val; % no scaling

        obj_val_eps_i_plus=compute_obj_val_eps_i_indirect_func(1,eps_val_scaled_i,param,x_cell,dx_dparam_i_cell,obj_func_given_x,other_input_obj_func0);

        obj_val_eps_i_minus=compute_obj_val_eps_i_indirect_func(-1,eps_val_scaled_i,param,x_cell,dx_dparam_i_cell,obj_func_given_x,other_input_obj_func0);

        dobj_dparam_indirect_i=(obj_val_eps_i_plus-obj_val_eps_i_minus)./(2*eps_val_scaled_i);

        if Richardson_spec==1
            %% Richardson extrapolation
            obj_val_eps_i_plus_half=compute_obj_val_eps_i_indirect_func(0.5,eps_val_scaled_i,param,x_cell,dx_dparam_i_cell,obj_func_given_x,other_input_obj_func0);

            obj_val_eps_i_minus_half=compute_obj_val_eps_i_indirect_func(-0.5,eps_val_scaled_i,param,x_cell,dx_dparam_i_cell,obj_func_given_x,other_input_obj_func0);

            A_h=(obj_val_eps_i_plus-obj_val_eps_i_minus)./(2*eps_val_scaled_i);
            A_h_half=(obj_val_eps_i_plus_half-obj_val_eps_i_minus_half)./(eps_val_scaled_i);

            dobj_dparam_indirect_i=(4*A_h_half-A_h)/3;
        end

        dobj_dparam_indirect(i,1)=dobj_dparam_indirect_i;
    end% for loop wrt i
end
