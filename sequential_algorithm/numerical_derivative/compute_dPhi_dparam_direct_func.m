%%%%%%%%%%%%%%%
function dPhi_dparam_direct_cell=...
    compute_dPhi_dparam_direct_func(x_cell,param,Phi_func,other_input_Phi_func0,Richardson_spec)


        n_param=size(param(:),1);
        eps_val=eps^(1/3); % central difference numerical derivative

    %% Compute dPhi_dparam_direct
        n_vars=size(x_cell,2);
    for m=1:n_vars
        dPhi_dparam_direct_cell{m}=[];
    end

    for m=1:n_vars
        size_x{m}=size(x_cell{m});
        cat_dim{m}=size(size_x{m},2)+1;
    end

    for i=1:n_param

        Phi_eps_i_plus_cell=compute_Phi_eps_i_direct_func(1,i,eps_val,param,x_cell,Phi_func,other_input_Phi_func0);

        Phi_eps_i_minus_cell=compute_Phi_eps_i_direct_func(-1,i,eps_val,param,x_cell,Phi_func,other_input_Phi_func0);

        for m=1:n_vars
            %%%dPhi_dparam_direct(:,i)=(Phi_eps_i_plus(:)-Phi_eps_i_minus(:))./(2*eps_val);
            dPhi_dparam_direct_m=(Phi_eps_i_plus_cell{m}-Phi_eps_i_minus_cell{m})./(2*eps_val);


            if Richardson_spec==1
                %% Richardson extrapolation
                Phi_eps_i_plus_half_cell=compute_Phi_eps_i_direct_func(0.5,i,eps_val,param,x_cell,Phi_func,other_input_Phi_func0);

                Phi_eps_i_minus_half_cell=compute_Phi_eps_i_direct_func(-0.5,i,eps_val,param,x_cell,Phi_func,other_input_Phi_func0);

                A_h=(Phi_eps_i_plus_cell{m}-Phi_eps_i_minus_cell{m})./(2*eps_val);
                A_h_half=(Phi_eps_i_plus_half_cell{m}-Phi_eps_i_minus_half_cell{m})./(eps_val);

                dPhi_dparam_direct_m=(4*A_h_half-A_h)/3;
            end
            dPhi_dparam_direct_cell{m}=cat(cat_dim{m},dPhi_dparam_direct_cell{m},dPhi_dparam_direct_m);
        end

    end% for loop wrt i

end
