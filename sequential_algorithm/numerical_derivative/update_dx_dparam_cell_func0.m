function [out,other_vars]=...
    update_dx_dparam_cell_func0(varargin)

        n_vars=varargin{end};

        for m=1:n_vars
            dx_dparam_initial_cell{m}=varargin{m};
            x_cell{m}=varargin{n_vars+m};
            Phi_cell{m}=varargin{2*n_vars+1+m};
        end

        param=varargin{2*n_vars+1};

        c=0;
        c=c+1;spec_gmres=varargin{end-c};
        c=c+1;solve_linear_eq_spec_dx_dparam=varargin{end-c};
        c=c+1;middle_vars_cell_plus=varargin{end-c};
        c=c+1;middle_vars_cell_minus=varargin{end-c};
        c=c+1;other_input_Phi_func=varargin{end-c};
        c=c+1;other_input_Phi_func0=varargin{end-c};
        c=c+1;Phi_func=varargin{end-c};
        c=c+1;central_difference_spec=varargin{end-c};
        c=c+1;Richardson_spec=varargin{end-c};
        c=c+1;eps_val=varargin{end-c};    
        
    for m=1:n_vars
        size_x{m}=size(x_cell{m});
        cat_dim{m}=size(size_x{m},2)+1;
    end


    if central_difference_spec==1
        eps_val=eps^(1/3);
    elseif central_difference_spec==0
        eps_val=sqrt(eps);
    end

    n_param=size(param(:),1);


    other_vars=[];

    %%%%%%%%%%%%%%%%%%%%
    %% Compute dPhi_dparam_direct
    dPhi_dparam_direct_cell=compute_dPhi_dparam_direct_func(x_cell,param,Phi_func,...
    other_input_Phi_func0,Richardson_spec);

    %%%%%%%%%%%%%%%
    if solve_linear_eq_spec_dx_dparam==0
        %% Compute dPhi_dparam_indirect
        for m=1:n_vars
            dx_dparam_initial_reshaped_cell{m}=reshape(dx_dparam_initial_cell{m},size(x_cell{m}(:),1),size(param(:),1));
        end
        
        norm_dx_dparam_i=zeros(n_vars,1);
        for m=1:n_vars
            dPhi_dparam_indirect_cell{m}=[];
        end

        for i=1:n_param
            
            for m=1:n_vars
                dx_dparam_i_initial_cell{m}=reshape(dx_dparam_initial_reshaped_cell{m}(:,i),size_x{m});
                norm_dx_dparam_i(m)=max(abs(dx_dparam_i_initial_cell{m}(:)));
            end
            eps_val_scaled_i=eps_val./max(max(norm_dx_dparam_i),1e-8);%scaling
            %eps_val_scaled_i=eps_val; % no scaling
            
            Phi_eps_i_plus_cell=compute_Phi_eps_i_indirect_func(1,eps_val_scaled_i,param,x_cell,dx_dparam_i_initial_cell,Phi_func,other_input_Phi_func0);

            Phi_eps_i_minus_cell=compute_Phi_eps_i_indirect_func(-1,eps_val_scaled_i,param,x_cell,dx_dparam_i_initial_cell,Phi_func,other_input_Phi_func0);

            for m=1:n_vars
                dPhi_dparam_indirect_m=(Phi_eps_i_plus_cell{m}-Phi_eps_i_minus_cell{m})./(2*eps_val_scaled_i);

                if Richardson_spec==1
                    %% Richardson extrapolation
                    Phi_eps_i_plus_half_cell=compute_Phi_eps_i_indirect_func(0.5,eps_val_scaled_i,param,x_cell,dx_dparam_i_initial_cell,Phi_func,other_input_Phi_func0);

                    Phi_eps_i_minus_half_cell=compute_Phi_eps_i_indirect_func(-0.5,eps_val_scaled_i,param,x_cell,dx_dparam_i_initial_cell,Phi_func,other_input_Phi_func0);

                    A_h=(Phi_eps_i_plus_cell{m}-Phi_eps_i_minus_cell{m})./(2*eps_val_scaled_i);
                    A_h_half=(Phi_eps_i_plus_half_cell{m}-Phi_eps_i_minus_half_cell{m})./(eps_val_scaled_i);

                    dPhi_dparam_indirect_m=(4*A_h_half-A_h)/3;
                end

                dPhi_dparam_indirect_cell{m}=cat(cat_dim{m},dPhi_dparam_indirect_cell{m},dPhi_dparam_indirect_m);
            end

        end% for loop wrt i


        for m=1:n_vars
            dx_dparam_updated_cell{m}=dPhi_dparam_direct_cell{m}+dPhi_dparam_indirect_cell{m};%n_x*n_param
        end

        
        %%resid_param=resid_param*0;%%%%
        %resid_dx_dparam=resid_dx_dparam*0;%%%%

    elseif solve_linear_eq_spec_dx_dparam==1
        %%% Solve for dx_dparam, by solving linear eq %%%
        fixed_point_iter_spec=1;

        [dx_dparam_updated_cell,flag,relres,iter,resvec]=...
                    solve_linear_eq_JF_cell_func(dPhi_dparam_direct_cell,dx_dparam_initial_cell,spec_gmres,Phi_func,n_vars,...
                   x_cell{:},param,other_input_Phi_func0{:},...
                        fixed_point_iter_spec);
        
        other_vars.iter_gmres=iter;
        other_vars.resvec=resvec;
        other_vars.flag=flag;
        
    end%solve_linear_eq_spec_dx_dparam==0 or 1
    
    out={dx_dparam_updated_cell{:}};

end
