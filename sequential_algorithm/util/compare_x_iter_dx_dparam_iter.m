%%% Input: 
% Phi_func, x_initial0_cell, param_initial,other_input_Phi_func

%% Compare x-iterations (nonlinear) and dx-dparam-iterations (linear; GMRES)
n_vars=size(x_initial0_cell,2);

%% Solve for x
spec_x_iter=[];
spec_x_iter.TOL=1e-12;
spec_x_iter.Anderson_acceleration=1;%%
[x_sol_cell,~,iter_info_x]=...
    spectral_func(Phi_func,spec_x_iter,{x_initial0_cell{:}},param_initial,other_input_Phi_func{:});
Phi_cell=Phi_func(x_sol_cell{:},param_initial,other_input_Phi_func{:});

eps_val=eps^(1/3);
if 1==1%% Perturb
    i=1;
    param_initial_eps_i=param_initial;
    param_initial_eps_i(i)=param_initial(i)+eps_val;
    [x_eps_sol_cell,~,iter_info_x_eps]=...
        spectral_func(Phi_func,spec_x_iter,{x_sol_cell{:}},param_initial_eps_i,other_input_Phi_func{:});
    
    param_initial_eps_i(i)=param_initial(i)-eps_val;
    [x_eps_minus_sol_cell,~,iter_info_x_eps_minus]=...
        spectral_func(Phi_func,spec_x_iter,{x_sol_cell{:}},param_initial_eps_i,other_input_Phi_func{:});
    
    dx_dparam_i_numer_diff=(x_eps_sol_cell{i}-x_eps_minus_sol_cell{i})/(2*eps_val);
    dx_dparam_numer_diff_vec_i=dx_dparam_i_numer_diff(:);
end
Phi_cell=Phi_func(x_sol_cell{:},param_initial,other_input_Phi_func{:});

%% Solve for dx_dparam
for ii=1:n_vars
    dx_dparam_initial0_cell{ii}=zeros([size(x_initial0_cell{ii}),size(param_initial(:),1)]);
end

central_difference_spec=1;
solve_linear_eq_spec_dx_dparam=1;

spec_gmres=[];
spec_gmres.ITER_MAX=200;
spec_gmres.TOL=1e-6;

Richardson_spec=0;

[dx_dparam_sol_cell,iter_gmres_dx_dparam]=...
    update_dx_dparam_cell_func0(dx_dparam_initial0_cell{:},...
    x_sol_cell{:},param_initial,Phi_cell{:},eps_val,...
        Richardson_spec,central_difference_spec,Phi_func,...
other_input_Phi_func,[],...
[],[],solve_linear_eq_spec_dx_dparam,spec_gmres,Richardson_spec,n_vars);

for ii=1:size(iter_gmres_dx_dparam.iter_gmres,2)
    iter_info_gmres_dx_dparam.feval(ii)=prod(iter_gmres_dx_dparam.iter_gmres{ii});
end
iter_info_gmres_dx_dparam.flag=iter_gmres_dx_dparam.flag;

dx_dparam_sol_vec=reshape(dx_dparam_sol_cell{1},[],size(param_initial(:),1));
dx_dparam_sol_vec_i=dx_dparam_sol_vec(:,i);

log10_comparison_mat=[log10(iter_gmres_dx_dparam.resvec{1}(1:50)),log10(iter_info_x_eps.DIST_table(1:50))];

