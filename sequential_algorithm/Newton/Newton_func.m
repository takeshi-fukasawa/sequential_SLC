function  [x_sol_cell,iter_info,other_vars,...
    flag_vec,relres_vec,iter_cell,resvec_cell]=...
    Newton_func(varargin)

quasi_newton_spec=0;

n_vars=varargin{end};
merit_func=varargin{end-1};
spec_gmres=varargin{end-2};

%%% varargin{n_vars+3:end-3}: other input for the function

for i=1:n_vars
   x_0_cell{i}=varargin{i};
   cell_size{i}=size(x_0_cell{i});
   elem_x(i)=prod(cell_size{i});
end
x_k_cell=x_0_cell;

n_x_all=sum(elem_x);

if 1==0
    B_k=eye(n_x_all);
end

fun=varargin{n_vars+1};%fixed point mapping
fixed_point_iter_spec=varargin{n_vars+2};


%%% Preliminary parameters %%%
sigma=0.8;
alpha=1e-6;
%alpha=1e-3;

if quasi_newton_spec==0
    ITER_MAX=40;
else
    ITER_MAX=200;
end

TOL=1e-10;
line_search_spec=1;%%%%
ITER_MAX_LINE_SEARCH=10;
flag_vec=[];
relres_vec=[];
iter_cell=[];
resvec_cell=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%

ITER_line_search_table=NaN(ITER_MAX,1);
step_size_table=NaN(ITER_MAX,1);
obj_val_table=NaN(ITER_MAX,1);
DIST_table=NaN(ITER_MAX,1);
feval=0;

%%%%%%%%%%%%%%
%spec_gmres.Newton_update_func=@compute_z_func_given_v_psteady_unobs_hetero;
%%%%%%%%%%%%%%%%

tic
for k=0:ITER_MAX-1

   if quasi_newton_spec==0% JFNK update 

    if isfield(spec_gmres,'Newton_update_func')==1
        [x_k_plus_1_cell,Delta_x_k_cell,F_k_vec,F_k_cell,other_vars]=spec_gmres.Newton_update_func(...
            x_k_cell{:},varargin{n_vars+3:end-3});

    else
        [x_k_plus_1_cell,Delta_x_k_cell,F_k_vec,F_k_cell,other_vars,...
            flag_vec,relres_vec,iter_cell,resvec_cell]=...
            Newton_update_func(x_k_cell{:},varargin{n_vars+1:end-3},spec_gmres,n_vars);
    end% compute_z_y_func

   if isempty(iter_cell)==0
     feval=feval+2*iter_cell{1}(2);
   else
       feval=feval+1;
   end

    if flag_vec~=0 & 1==0
        message="GMRES not converge?"
        [k,flag_vec,iter_cell]
        resvec_cell{1}(end)
    end

    else % Broyden update

           [out_temp,other_vars]=fun(x_k_cell{:},...
                varargin{n_vars+3:end-3});
           feval=feval+1;

           if fixed_point_iter_spec==1
               for i=1:n_vars
                   F_k_cell{i}=out_temp{i}-x_k_cell{i};
               end
           else%fixed_point_iter_spec==0
               x_k_vec=cell_to_vec_func(x_k_cell);
               F_k_cell=out_temp;
           end
           F_k_vec=cell_to_vec_func(F_k_cell);

          Delta_x_k_vec=-(inv(B_k))*F_k_vec;
          Delta_x_k_cell=vec_to_cell_func1(Delta_x_k_vec,cell_size);

     end% Broyden update or not

    merit_val_k=merit_func(F_k_cell);
    obj_val_table(k+1,1)=merit_val_k;

    if line_search_spec==1
        step_size=1;
        for iter_line_search=1:ITER_MAX_LINE_SEARCH
           for i=1:n_vars
              x_k_plus_1_cell{i}=x_k_cell{i}+step_size*Delta_x_k_cell{i};
           end
    
           [out_temp,other_vars]=fun(x_k_plus_1_cell{:},...
                varargin{n_vars+3:end-3});

           feval=feval+1;

           if fixed_point_iter_spec==1
               for i=1:n_vars
                   F_k_plus_1_cell{i}=out_temp{i}-x_k_plus_1_cell{i};
               end
           else%fixed_point_iter_spec==0
               F_k_plus_1_cell=out_temp;
           end
     
            merit_val_k_plus_1=merit_func(F_k_plus_1_cell);
    
            LHS=merit_val_k_plus_1;

            if quasi_newton_spec==0
                eta_k=1/((k+1)^2);
                RHS=(1-1*alpha*step_size)*merit_val_k+eta_k*0.00;
            else% quasi_newton_dummy==1
                eta_k=1/((k+1)^2);
                RHS=(1+eta_k)*merit_val_k-alpha*(step_size^2)*sum(Delta_x_k_vec.^2);
            end %quasi_Newton_dummy==0 or 1

            %[LHS,RHS]
            terminate_dummy=(LHS<=RHS);
            %terminate_dummy=1;%%%%%%

            if terminate_dummy==1
                break;
            else
                step_size=step_size*sigma;
            end% Assess the termination of line search
        end%  for iter_line_search=1:ITER_MAX_LINE_SEARCH
        
        ITER_line_search_table(k+1,1)=iter_line_search;
        step_size_table(k+1,1)=step_size;

    end% line_search_spec==1
    DIST=sup_norm_func(F_k_cell);
    DIST_table(k+1,1)=DIST;
    %[merit_val_k_plus_1,merit_val_k]
    
    if DIST>TOL
       if quasi_newton_spec==1
            s_k_vec=cell_to_vec_func(x_k_plus_1_cell)-cell_to_vec_func(x_k_cell);
            y_k_vec=cell_to_vec_func(F_k_plus_1_cell)-cell_to_vec_func(F_k_cell);
            B_k=B_k+(y_k_vec-B_k*s_k_vec)*(s_k_vec')./sum(s_k_vec.^2);
       end % quasi_newton_spec==1

       x_k_cell=x_k_plus_1_cell;

    else
        break;
    end

end% k=0:ITER_MAX
iter_info.t_cpu=toc;

x_sol_cell=x_k_plus_1_cell;

iter_info.obj_val_table=obj_val_table;
iter_info.DIST_table=DIST_table;%L2 norm
iter_info.ITER_line_search_table=ITER_line_search_table;
iter_info.step_size_table=step_size_table;
iter_info.n_iter=k+1;
iter_info.feval=feval;
iter_info.obj_val=merit_val_k;


end



