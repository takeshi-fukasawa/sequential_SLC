function obj_val_eps_i=compute_obj_val_eps_i_indirect_func(val,eps_val_scaled_i,param,x_cell,dx_dparam_i_cell,obj_func_given_x,other_input_obj_func0)

    n_vars=size(x_cell,2);
    for m=1:n_vars
        x_eps_i_cell{m}=x_cell{m}+val*eps_val_scaled_i*dx_dparam_i_cell{m};
    end

    if isempty(other_input_obj_func0)==1
        obj_val_eps_i=obj_func_given_x(param,x_eps_i_cell{:});
    else
        obj_val_eps_i=obj_func_given_x(param,x_eps_i_cell{:},other_input_obj_func0{:});
    end
end
