function [J,cell_size,n_elem_temp]=...
    compute_Jacobian_func(X_cell,cell_func,fixed_point_iter_spec,varargin)


    %rng(100)

        X_vec=cell_to_vec_func(X_cell);
        %%X_vec=X_vec.*exp(1.5*rand(size(X_vec)));

        n=size(X_vec,1);

        n_vars_all=size(X_cell,2);

        n_elem_temp=zeros(n_vars_all,1);
        for i=1:n_vars_all
            cell_size{i}=size(X_cell{i});
            n_elem_temp(i)=prod(cell_size{i});
        end

        base_vec=vec_func(X_vec,cell_func,cell_size,...
            varargin{:});
            
        
        J=[];
        perturb_val=eps^(1/3);
        for i=1:n
            X_vec_plus_eps=X_vec;
            X_vec_plus_eps(i)=X_vec(i)+perturb_val;
            perturb_vec_plus_i=vec_func(X_vec_plus_eps,cell_func,cell_size,...
            varargin{:});

            X_vec_minus_eps=X_vec;
            X_vec_minus_eps(i)=X_vec(i)-perturb_val;
            perturb_vec_minus_i=vec_func(X_vec_minus_eps,cell_func,cell_size,...
            varargin{:});

            %%J=[J,(perturb_vec_plus_i-base_vec)./perturb_val];
            J=[J,(perturb_vec_plus_i-perturb_vec_minus_i)./(2*perturb_val)];
        end

        if fixed_point_iter_spec==1
            J=eye(n)-J;
        end

end
