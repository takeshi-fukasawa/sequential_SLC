function dobj_dparam=compute_dobj_dparam_func(param_initial,x_initial_cell,...
    eps_val,central_difference_spec,Phi_func,other_input_Phi_func0,other_input_Phi_func,...
            spec_gmres,...
        obj_func_given_x,other_input_obj_func0,Richardson_spec)

        n_vars_x=size(x_initial_cell,2);
        for m=1:n_vars_x
            size_x{m}=size(x_initial_cell{m});
        end
    
        
    dobj_dparam_direct=compute_dobj_dparam_direct_func(...
        param_initial,x_initial_cell,obj_func_given_x,other_input_obj_func0,Richardson_spec);

    if isempty(other_input_Phi_func)==1
        out=Phi_func(x_initial_cell{:},param_initial);
     else
        out=Phi_func(x_initial_cell{:},param_initial,other_input_Phi_func{:});
     end
 
     Phi_cell=out;
 
    partially_solve_spec=1;%%%
    [dx_dparam_cell,~]=update_dx_dparam_cell_func0([],x_initial_cell{:},param_initial,...
        Phi_cell,eps_val,Richardson_spec,central_difference_spec,Phi_func,other_input_Phi_func0,other_input_Phi_func,...
            [],[],...
            partially_solve_spec,spec_gmres,n_vars_x);

    dobj_dparam_indirect=compute_dobj_dparam_indirect_func(...
        param_initial,x_initial_cell,dx_dparam_cell,size_x,obj_func_given_x,other_input_obj_func0,Richardson_spec,n_vars_x);

    dobj_dparam=dobj_dparam_direct+dobj_dparam_indirect;

end
