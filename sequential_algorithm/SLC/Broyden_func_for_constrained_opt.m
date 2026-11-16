function [x_sol_cell,other_output_k_plus_1,iter_info]=...
    Broyden_func_for_constrained_opt(fun,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa.
global k_iter

n_var=size(x_0_cell,2);
mu_val_k=0.1;

spec=preliminary_setting_func(spec,n_var);

other_input_cell=varargin;

for i=1:n_var
    cell_size{i}=size(x_0_cell{1,i});
end

FLAG_ERROR=0;
k_iter=0;

x_max_cell=spec.x_max_cell;
x_min_cell=spec.x_min_cell;

obj_func=spec.obj_func;
Phi_func=spec.Phi_func;
other_input_obj_func=spec.other_input_obj_func;
other_input_Phi_func=spec.other_input_Phi_func;
param_max=spec.param_max;
param_min=spec.param_min;

% varargin:1*XXX

%% Read inputs
ITER_MAX=spec.ITER_MAX;
DIST_table=NaN(ITER_MAX,n_var);
ITER_table_LINE_SEARCH=NaN(ITER_MAX,1);
step_size_table=NaN(ITER_MAX,n_var);
obj_val_table=NaN(ITER_MAX,n_var);
merit_val_table=NaN(ITER_MAX,3);
mu_val_table=NaN(ITER_MAX,1);

tic_spectral=tic;
 if spec.bound_spec==1
    x_0_cell=projection_func(x_0_cell,x_max_cell,x_min_cell);
 end
 [fun_0_cell_temp,other_output_0]=fun(x_0_cell{:},other_input_cell{:});

if spec.fixed_point_iter_spec==1
     for i=1:n_var
         fun_0_cell{1,i}=fun_0_cell_temp{1,i}-x_0_cell{1,i};
     end
else
       fun_0_cell=fun_0_cell_temp;
end


feval=1;
    
    %%% DIST: norm of F(x)=x-Phi(x). 
    DIST_vec=ones(1,n_var);
    for i=1:n_var
      DIST_vec(1,i)=norm_func(fun_0_cell{1,i}(:),x_0_cell{1,i}(:),spec.norm_spec(i));
      obj_val_table(1,i)=sum(fun_0_cell{1,i}(:).^2);%(L2 norm)^2

      spectral_coef_0{1,i}=spec.spectral_coef_0;
      
      if spec.update_spec==0
        spectral_coef_0{1,i}=1;
      end

      spectral_coef_table(1,i)=spectral_coef_0{1,i};
      
    end % loop wrt i


    DIST=nanmax(DIST_vec);
    DIST_table(1,:)=DIST_vec;

    merit_val_vec_0=compute_merit_val_func(x_0_cell,mu_val_k,...
        obj_func,Phi_func,...
        other_input_obj_func,other_input_Phi_func,...
            spec.fixed_point_iter_spec,param_max,param_min);


conv=(sum((DIST_vec<spec.TOL),'all')==n_var);

x_k_cell=x_0_cell;
fun_k_cell=fun_0_cell;
other_output_k_plus_1=other_output_0;
merit_val_vec_k=merit_val_vec_0;

x_k_minus_1_cell=x_0_cell;%%%
fun_k_minus_1_cell=fun_0_cell;%%%

s_mat=[];
y_mat=[];

%%%%%%%% Loop %%%%%%%%%%%

if conv==0 & ITER_MAX>=2

for k=0:ITER_MAX-2
    k_iter=k;
    
    if k>=1

        s_vec=[];
        y_vec=[];
        for i=1:n_var
            Delta_x_cell{1,i}=x_k_cell{i}-x_k_minus_1_cell{i};
            Delta_fun_cell{1,i}=fun_k_cell{i}-fun_k_minus_1_cell{i};
            s_vec=[s_vec;Delta_x_cell{1,i}(:)];
            y_vec=[y_vec;Delta_fun_cell{1,i}(:)];
        end % loop wrt i

        s_mat=[s_mat,s_vec];
        y_mat=[y_mat,y_vec];
        
        dummy=(isempty(spec.update_spec)==1 || spec.update_spec~=0);
        
        if dummy==1
            [spectral_coef_k,spectral_coef_max]=compute_spectral_coef_func(...
            Delta_x_cell,Delta_fun_cell,spec,k);
            spec.spectral_coef_max=spectral_coef_max;
        end

    else% k==0
        for i=1:n_var
            spectral_coef_k{1,i}=spectral_coef_0{1,i};
            if isempty(spec.dampening_param)==0
                spectral_coef_k{1,i}=spectral_coef_k{1,i}*(spec.dampening_param{1,i});
            end
        end
    end% k>=1 or 0

   
    if spec.update_spec==0 %| k_iter>=30
        for i=1:n_var
            spectral_coef_k{1,i}=1;
        end
    end  %if spec.update_spec==0


    for i=1:n_var
        spectral_coef_vec(1,i)=max(spectral_coef_k{1,i}(:));
    end


    %%% Update variables %%%%%%%%%%%%%%%

    for i=1:n_var
        d_k_cell{1,i}=fun_k_cell{1,i};
    end % for loop wrt i

    %%%% Broyden update %%%%
    if 1==1
        for i=1:n_var
            d_k_cell_spectral{1,i}=d_k_cell{1,i}.*spectral_coef_k{1,i};
        end
        
        z=cell_to_vec_func(fun_k_cell);

        if k==0
            s_mat=[-cell_to_vec_func(fun_k_cell)];%n_elem*1
        end

        small_val=1e-4;

        if k>=1
            j_start=1;
            j_end=k;
            if k>=20 & 1==1
                % Avoid using too old (sometimes unreliable) information
                j_start=20;%%%%
                j_end=k;%%%%

                %j_start=1;
                %j_end=20;
            end

            for j=j_start:j_end
                denom=sum(s_mat(:,j).^2,'omitnan');
                denom=max(small_val,denom);% Avoid division by too small number
                temp=sum(s_mat(:,j).*z,'omitnan');%scalar
                z=z+(s_mat(:,j+1))*temp./denom;
            end
        end

        denom=sum(s_mat(:,k+1).^2,'omitnan');
        denom=max(small_val,denom);% Avoid division by too small number

        temp=sum(s_mat(:,k+1).*z,'omitnan');%scalar
        s_k_plus_1=z./(1-temp./denom);

        d_k_cell=vec_to_cell_func1(s_k_plus_1,cell_size);

        for i=1:n_var
            spectral_coef_k{1,i}=1;
        end

        if k<=5
            d_k_cell=d_k_cell_spectral;
        end
        
    end%1==0
    %%%%%%%%%%%%%%%

    [x_k_plus_1_cell, fun_k_plus_1_cell,...
    other_output_k_plus_1,DIST_vec,iter_line_search,...
    obj_val_vec,step_size,mu_val_k_plus_1,merit_val_vec_k_plus_1]=...
        update_func_Broyden_for_constrained_opt(fun,x_k_cell,d_k_cell,d_k_cell_spectral,spectral_coef_k,mu_val_k,merit_val_vec_k,other_input_cell,...
        n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);

    feval=feval+iter_line_search;

    ITER_table_LINE_SEARCH(k+1,1)=iter_line_search;%% Number of line search iterations
    DIST_table(k+2,:)=DIST_vec;
    spectral_coef_table(k+1,:)=spectral_coef_vec;
    obj_val_table(k+2,:)=obj_val_vec;
    merit_val_table(k+2,:)=merit_val_vec_k_plus_1;
    mu_val_table(k+2,:)=mu_val_k_plus_1;
    step_size_table(k+1,:)=step_size;
    DIST=max(DIST_vec);
            

    if isnan(sum(DIST_vec))==1||isinf(sum(DIST_vec))==1
       message="NaN or Inf DIST"
       x_k_plus_1_cell=x_k_cell;
       FLAG_ERROR=1;
       break;
    end
   
   interval=10;
    if k-floor(k/interval)*interval==0&spec.DEBUG==1
        DIST_vec
    end

    if sum((DIST_vec<spec.TOL),'all')==n_var
        FLAG_ERROR=0;
        %DIST
        break;
    end

    %%% Replace variables for the next iteration
    x_k_minus_1_cell=x_k_cell;
    fun_k_minus_1_cell=fun_k_cell;


	x_k_cell=x_k_plus_1_cell;
    fun_k_cell=fun_k_plus_1_cell;

    merit_val_vec_k=merit_val_vec_k_plus_1;
    mu_val_k=mu_val_k_plus_1;

end %% end of for loop wrt k=0:ITER_MAX-1


else % no iteration
    k=0;
    if spec.fixed_point_iter_spec==1
        x_k_plus_1_cell=fun_0_cell_temp;% Use f(x0), where f is the fixed point mapping
    else
       x_k_plus_1_cell=x_0_cell;% Use x0
    end

    other_output_k_plus_1=other_output_0;
    fun_k_cell=fun_0_cell;
end

%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc(tic_spectral);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=feval;
iter_info.ITER_MAX=spec.ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;
iter_info.obj_val_table=obj_val_table;
iter_info.merit_val_table=merit_val_table;
iter_info.mu_val_table=mu_val_table;

iter_info.fun_cell=fun_k_cell;

iter_info.DIST_table=DIST_table;
iter_info.spectral_coef_table=spectral_coef_table;
iter_info.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;

iter_info.step_size_table=step_size_table;
iter_info.spec=spec;


iter_info.fun=fun;

end

