function y = JV_APPROX_cell_func(v, F,n_vars,varargin)

    %% Compute Jacovian-vector product (dG/dx)*v, where G is the function for which we want to solve G(x)=0. This function is used in the GMRES algorithm to solve linear equations involving the Jacobian matrix dG/dx without explicitly computing the Jacobian matrix itself. See solve_linear_eq_JF_cell_func.m for details.
    
    %%% NaN values: Ignored in the GMRES algorithm. 

    % v: vector
    % F: function handle for G

    for i=1:n_vars
        x_cell{i}=varargin{i};% cell
        cell_size{i}=size(x_cell{i});
    end

    
    fixed_point_iter_spec=varargin{end-3};
    Fx_cell=varargin{end-2};
    scaling_cell=varargin{end-1};
    nonnan_ids=varargin{end};

    n_elem=0;
    for i=1:n_vars
        n_elem=n_elem+prod(cell_size{i});
    end
    v_temp=NaN(n_elem,1);
    v_temp(nonnan_ids)=v;
    v_cell=vec_to_cell_func1(v_temp,cell_size);
 
    central_difference_spec=1;

    if central_difference_spec==0
        eps_val=sqrt(eps);
    else
        eps_val=(eps.^(1/3));
    end
    eps_val=eps_val./max(max(abs(v(:))),1e-8);
    eps_val=min(eps_val,1e-4);%% Avoid too large eps_val
    
    if isempty(Fx_cell)==1 & central_difference_spec==0
        Fx_cell=G_func(F,x_cell,fixed_point_iter_spec,varargin{1,1+n_vars:end-4});% unperturbed residual

        for m=1:n_vars
            Fx_per_cell{m}=Fx_per_cell{m}./scaling_cell{m};
        end
    end
    
    x_eps_plus_cell=x_cell;
     for i=1:n_vars
           x_eps_plus_cell{i}=x_cell{i}+eps_val*v_cell{i};
      end

    [Fx_per_plus_cell]=G_func(F,x_eps_plus_cell,fixed_point_iter_spec,varargin{1,1+n_vars:end-4});%perturbed residual

    for m=1:n_vars
        Fx_per_plus_cell{m}=Fx_per_plus_cell{m}./scaling_cell{m};
    end

    Fx_per_plus_vec=cell_to_vec_func(Fx_per_plus_cell);

    if central_difference_spec==0
        y = (Fx_per_plus_vec - cell_to_vec_func(Fx_cell)) / eps_val; % approximation of jacobian action on krylov vector
    else

        x_eps_minus_cell=x_cell;
        for i=1:n_vars
           x_eps_minus_cell{i}=x_cell{i}-eps_val*v_cell{i};
        end

        [Fx_per_minus_cell]=G_func(F,x_eps_minus_cell,fixed_point_iter_spec,varargin{1,1+n_vars:end-4});

        for m=1:n_vars
            Fx_per_minus_cell{m}=Fx_per_minus_cell{m}./scaling_cell{m};
        end

        y = (Fx_per_plus_vec - cell_to_vec_func(Fx_per_minus_cell)) / (2*eps_val); 
    end

    y=y(isnan(y)==0);% Omit NaN values
    y=y(:);% vector

end