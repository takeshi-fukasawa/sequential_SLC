function Phi_eps_i_cell=compute_Phi_eps_i_direct_func(val,i,eps_val,param,x_cell,Phi_func,other_input_Phi_func0)

    param_eps_i=param;       
    param_eps_i(i)=param(i)+val*eps_val;

    if isempty(other_input_Phi_func0)==1       
        out=Phi_func(x_cell{:},param_eps_i);
    else
        out=Phi_func(x_cell{:},param_eps_i,other_input_Phi_func0{:});
    end
    Phi_eps_i_cell=out;
end