function [x_k_plus_1_cell,Delta_x_k_cell,Fx_k,F_k_cell,other_vars,...
    flag_vec,relres_vec,iter_cell,resvec_cell]=Newton_update_func(varargin)

    %%% Use the idea of Jacobian-free Krylov
    %%% Idea: Solve (F'(x))d=b for d
     n_vars=varargin{end};
     spec_gmres=varargin{end-1};
     for i=1:n_vars
        x_k_cell{i}=varargin{i};
     end
    fun=varargin{n_vars+1};%fixed point mapping
    fixed_point_iter_spec=varargin{n_vars+2};

    F_k_cell=G_func(fun,x_k_cell,fixed_point_iter_spec,varargin{n_vars+3:end-2});

    Fx_k=cell_to_vec_func(F_k_cell);

    
    if 1==1
        v_init=[];

        [d_cell,flag_vec,relres_vec,iter_cell,resvec_cell]=solve_linear_eq_JF_cell_func(...
            F_k_cell,v_init,spec_gmres,fun,n_vars,...
        x_k_cell{:},varargin{n_vars+3:end-2},fixed_point_iter_spec);
        
        d=[];
        for i=1:n_vars
            d=[d;d_cell{i}(:)];
        end

        %if flag_vec~=0
        %    message="Nonconvergence of GMRES"
        %end
    else
        %%%% Compute Jacobian (temp) %%%
        [J,cell_size,n_elem_temp]=...
            compute_Jacobian_func(x_k_cell,fun,fixed_point_iter_spec,varargin{n_vars+3:end-2});
    
        %v=Fx_k;
        
        %JV_approx=j_vec_approx(v);
        
        %JV=J*v;
    
        %difference_JV=JV_approx-JV;

        %%% Using Jacobian
        v_init=[];
        b=Fx_k;
        [sol,flag,relres,iter,resvec] = gmres(J,b,[],1e-7,100,[],[],v_init);
        %flag
        %iter
        %resid_norm=norm(J*sol-b);
        
        %resid_norm_JF=norm(J*d-Fx_k)

        flag_vec=[];
        relres_vec=[];
        iter_cell=[];
        resvec_cell=[];
        d=sol;%%%
    %%%%%%%%%%%%%%%%%%%%%
    end


    id=1;
    for i=1:n_vars
        d_i=d(id:id+numel(x_k_cell{i})-1);
         Delta_x_k_cell{i}=-reshape(d_i,size(x_k_cell{i}));        
         x_k_plus_1_cell{i}=x_k_cell{i}-reshape(d_i,size(x_k_cell{i}));
        id=id+numel(x_k_cell{i});
    end
 

    other_vars=[];
end
