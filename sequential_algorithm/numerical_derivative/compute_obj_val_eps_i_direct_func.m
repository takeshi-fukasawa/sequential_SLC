function obj_eps_i=compute_obj_val_eps_i_direct_func(val,i,eps_val,param,x_cell,obj_func_given_x,other_input_obj_func0)

    param_eps_i=param;       
    param_eps_i(i)=param(i)+val*eps_val;

    if isempty(other_input_obj_func0)==1       
        out=obj_func_given_x(param_eps_i,x_cell{:});
    else
        out=obj_func_given_x(param_eps_i,x_cell{:},other_input_obj_func0{:});
    end
    obj_eps_i=out;
end