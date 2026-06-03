function [G_cell]=G_func(fun,x_cell,fixed_point_iter_spec,varargin)

    %% Compute G(x_cell), where G is the function for which we want to solve G(x)=0. If fixed_point_iter_spec==1, G(x)=x-fun(x). If fixed_point_iter_spec==0, G(x)=fun(x). 

    
    global G_scaling_cell


    n_vars=size(x_cell,2);
    G_cell=fun(x_cell{:},varargin{:});

    if fixed_point_iter_spec==1
        for i=1:n_vars
            G_cell{i}=x_cell{i}-G_cell{i};
        end
    end

    %%% Scaled, in order to avoid ill-conditioned GMRES
    for i=1:n_vars
        G_cell{i}=G_cell{i}./G_scaling_cell{i};
    end

end
