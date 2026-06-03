function [v_sol_cell,flag_vec,relres_vec,iter_cell,resvec_cell]=...
    solve_linear_eq_JF_cell_func(b_cell,v_cell_init,spec_gmres,varargin)

    %%% Solve (dG/dx)^(-1) vec=b for vec.
    %%% Jacobian-free method.
    %%% If fixed_point_iter_spec==1, solve (I - dPhi/dx)^(-1) vec=b.
    
    %fun=varargin{1};
    %n_vars=varargin{2};
    %x=varargin{3};
    
    K=numel(b_cell{1})/numel(varargin{3});%Number of columns in b_cell{m} (number of different b's to solve for)
    
     if isfield(spec_gmres,'fixed_point_iter_spec')==1
        fixed_point_iter_spec=spec_gmres.fixed_point_iter_spec;
     else
         fixed_point_iter_spec=0;
     end
    %fixed_point_iter_spec=varargin{end};
    
    n_vars=size(b_cell,2);

    b=[];
    for m=1:n_vars
        temp=std(b_cell{m}(:),'omitnan');
        if temp==0
            temp=1;
        end
       scaling_cell{m}=temp;

       scaling_cell{m}=1;%%%
       
       b_cell{m}=b_cell{m}./scaling_cell{m};

       b=[b;reshape(b_cell{m},[],K)];%[]*K
       cell_size{m}=[numel(b_cell{m})/K,1];
    end
    
    nonnan_ids=find(isnan(b(:,1))==0 & isinf(b(:,1))==0);

    %scaling_cell

    j_vec_approx=@(vec)JV_APPROX_cell_func(vec,varargin{:},[],...
        scaling_cell,nonnan_ids);

    v_init=[];
    if isempty(v_cell_init)==0
       for m=1:n_vars
          v_init=[v_init;reshape(v_cell_init{m},[],K)];%[]*K
       end
    end

    vec_sol=b*0;
    
    flag_vec=zeros(1,K);
    relres_vec=zeros(1,K);
    iter_cell=[];
    resvec_cell=[];

    if isfield(spec_gmres,'TOL')==0
        TOL=1e-8;
    else
        TOL=spec_gmres.TOL;
    end

    if isfield(spec_gmres,'ITER_MAX')==0
        ITER_MAX=100;
    else
        ITER_MAX=spec_gmres.ITER_MAX;
    end
    
    for i=1:K
        b_i=b(:,i);
        b_i=b_i(nonnan_ids);% Omit nan values
        if norm(b_i)<TOL
            vec_sol(:,i)=0;
            flag=[];relres=[];iter=[];resvec=[];
        else
            if isempty(v_init)==0
                v_init_i=v_init(nonnan_ids,i);
            else
                v_init_i=[];
            end
            restart=[];

            try
                 [vec_sol_i,flag_vec(i),relres_vec(i),iter_cell{i},resvec_cell{i}]=...
                gmres(...
                j_vec_approx,b_i,restart,TOL/norm(b_i),ITER_MAX,[],[],v_init_i);
                vec_sol(nonnan_ids,i)=vec_sol_i;
            catch
                vec_sol_i=b_i*0;
                vec_sol(nonnan_ids,i)=vec_sol_i;
                flag_vec(i)=NaN;relres_vec(i)=NaN;iter_cell{i}=NaN;resvec_cell{i}=NaN;

            end
        end%norm(b_i)<TOL or not
    end% for i=1:K

    
    for m=1:n_vars
        v_sol_cell{m}=0*zeros([cell_size{m}(1),K]);
    end
    for i=1:K
        v_sol_i_cell=vec_to_cell_func1(vec_sol(:,i),cell_size);
        for m=1:n_vars
            v_sol_cell{m}(:,i)=v_sol_i_cell{m}(:);%each cell is an []*K-dim matrix
        end%m=1:n_vars
    end%i=1:K

    for m=1:n_vars
        v_sol_cell{m}=reshape(v_sol_cell{m},size(b_cell{m}));
    end

end
