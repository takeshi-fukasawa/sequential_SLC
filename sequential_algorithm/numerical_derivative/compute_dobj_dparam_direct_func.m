%%%%%%%%%%%%%%%
function dobj_dparam_direct=compute_dobj_dparam_direct_func(param,x_cell,obj_fun,other_input_obj_func,Richardson_spec)

        dobj_dparam_direct=zeros(size(param));
        n_param=size(param(:),1);
        eps_val=eps^(1/3); % central difference numerical derivative
        for i=1:n_param

            obj_eps_i_plus=compute_obj_val_eps_i_direct_func(1,i,eps_val,param,x_cell,obj_fun,other_input_obj_func);

            obj_eps_i_minus=compute_obj_val_eps_i_direct_func(-1,i,eps_val,param,x_cell,obj_fun,other_input_obj_func);

            dobj_dparam_direct(i)=(obj_eps_i_plus-obj_eps_i_minus)./(2*eps_val);

            if Richardson_spec==1
                %% Richardson extrapolation 
                obj_eps_i_plus_half=compute_obj_val_eps_i_direct_func(0.5,i,eps_val,param,x_cell,obj_fun,other_input_obj_func);

                obj_eps_i_minus_half=compute_obj_val_eps_i_direct_func(-0.5,i,eps_val,param,x_cell,obj_fun,other_input_obj_func);
    
    
                A_h=(obj_eps_i_plus-obj_eps_i_minus)./(2*eps_val);
                A_h_half=(obj_eps_i_plus_half-obj_eps_i_minus_half)./(eps_val);
    
                dobj_dparam_direct(i)=(4*A_h_half-A_h)/3;
            end
        end% for loop wrt i
end
