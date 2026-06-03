function [ll,s0]=obj_nfxp_func_bhhh(...
    param,Phi_func,obj_func,...
other_vars_Phi_func,other_vars_obj_func,spec,varargin)

global x_initial_cell_global count_Phi_eval
global x_sol_cell
param'

x_initial0_cell=varargin;

if isempty(x_initial_cell_global)==0
    x_initial_cell=x_initial_cell_global;%hot start
else
    x_initial_cell=x_initial0_cell;
end

%%x_initial_cell=x_initial0_cell;

if isempty(spec)==1
    spec=[];
    spec.TOL=1e-12;
    spec.ITER_MAX=2000;
    spec.Anderson_acceleration=1;
end

spec.TOL=1e-12;%%%%%
spec.ITER_MAX=1000;%%%%%
spec.Anderson_acceleration=1;%%%%

[output_spectral,other_vars,iter_info]=...
        spectral_func(Phi_func,spec,...
        x_initial_cell,param,other_vars_Phi_func{:});

    if iter_info.FLAG_ERROR==1 | iter_info.feval>=iter_info.ITER_MAX
        temp="Not converge"
        iter_info.feval
        iter_info.DIST_table(end,:)
    end

    x_sol_cell=output_spectral;
    x_initial_cell_global=output_spectral;

    %iter_info.feval
    count_Phi_eval=count_Phi_eval+iter_info.feval;

    ll_i=obj_func(param,x_sol_cell{:},other_vars_obj_func{:});
    
    ll=mean(ll_i);
    
    %if spec.analytical_spec==0
    
    %% Compute gradient by numerical differentiation
    
        n_param=size(param,1);
        eps_val=eps^(1/3);
        n_obs=size(ll_i(:),1);
        s0=zeros(n_obs,n_param);

        for i=1:n_param
            param_eps_i=param;
            param_eps_i(i)=param(i)+eps_val;

            [output_spectral_eps_i,other_vars_eps_i,iter_info_eps_i]=...
            spectral_func(Phi_func,spec,...
            x_sol_cell,param_eps_i,other_vars_Phi_func{:});
            ll_i_eps=obj_func(param_eps_i,output_spectral_eps_i{:},other_vars_obj_func{:});

            param_eps_i_minus=param;
            param_eps_i_minus(i)=param(i)-eps_val;

            [output_spectral_eps_i,other_vars_eps_i,iter_info_eps_i]=...
            spectral_func(Phi_func,spec,...
            x_sol_cell,param_eps_i_minus,other_vars_Phi_func{:});
            ll_i_eps_minus=obj_func(param_eps_i_minus,output_spectral_eps_i{:},other_vars_obj_func{:});

            s0(:,i)=reshape(ll_i_eps-ll_i_eps_minus,n_obs,1)./(2*eps_val);

        end

    %else% analytical derivatives

    %end

    
end
