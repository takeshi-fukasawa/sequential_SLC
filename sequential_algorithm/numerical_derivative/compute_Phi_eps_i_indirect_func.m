function Phi_eps_i_cell=compute_Phi_eps_i_indirect_func(val,eps_val_scaled_i,param,x_cell,dx_dparam_i_cell,Phi_func,other_input_Phi_func0)

    n_vars=size(x_cell,2);
    for m=1:n_vars
        x_eps_i_cell{m}=x_cell{m}+val*eps_val_scaled_i*dx_dparam_i_cell{m};
    end

    if isempty(other_input_Phi_func0)==1
        out=Phi_func(x_eps_i_cell{:},param);
    else
        out=Phi_func(x_eps_i_cell{:},param,other_input_Phi_func0{:});
    end

    Phi_eps_i_cell=out;

end
