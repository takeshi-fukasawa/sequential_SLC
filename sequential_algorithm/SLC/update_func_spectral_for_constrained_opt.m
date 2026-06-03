%%%%%%%%%%%%%%%%%%
function [x_k_plus_1_cell, fun_k_plus_1_cell,...
other_output_k_plus_1,DIST_vec,iter_line_search,...
obj_val_vec,step_size,mu_val_k_plus_1,merit_val_vec_k_plus_1]=...
update_func_spectral_for_constrained_opt(update_fun,x_k_cell,d_k_cell,spectral_coef_k,mu_val_k,merit_val_vec_k,other_input_cell,...
n_var,spec,...
x_max_cell,x_min_cell,k,obj_val_table)

    step_size=ones(1,n_var);
    
    obj_func=spec.obj_func;
    Phi_func=spec.Phi_func;
    other_input_obj_func=spec.other_input_obj_func;
    other_input_Phi_func=spec.other_input_Phi_func;
    param_min=spec.param_min;
    param_max=spec.param_max;

    %%% Update variables
    if spec.line_search_spec==1
       spec.ITER_MAX_LINE_SEARCH=5;%%%
    else
        spec.ITER_MAX_LINE_SEARCH=1;
    end

    for iter_line_search=1:spec.ITER_MAX_LINE_SEARCH
       for i=1:n_var
            x_k_plus_1_cell{1,i}=x_k_cell{1,i}+step_size(i)*spectral_coef_k{1,i}.*d_k_cell{1,i};
       end % for loop wrt i

       if spec.bound_spec==1
           x_k_plus_1_cell=projection_func(x_k_plus_1_cell,x_max_cell,x_min_cell);

          % Adjust spctral_coef_k so that x_k_plus_1 is inside the boundary
          spectral_coef_k_temp=(x_k_plus_1_cell{1}-x_k_cell{1})./d_k_cell{1};%vector
           for i=1:n_var
              spectral_coef_k{i}=nanmax(spectral_coef_k_temp);
           end
       end

       if iter_line_search==1
            [fun_k_plus_1_cell,other_output_k_plus_1]=...
            update_fun(x_k_plus_1_cell{:},other_input_cell{:});
        

            if spec.fixed_point_iter_spec==1
                for i=1:n_var
                    fun_k_plus_1_cell{1,i}=fun_k_plus_1_cell{1,i}-x_k_plus_1_cell{1,i};
                end
            end
        end%iter_line_search==1

        %%% DIST: sup norm of F(x)=x-Phi(x). 
        DIST_vec=ones(1,n_var);
        obj_val_vec=zeros(1,n_var);

        for i=1:n_var
            DIST_vec(1,i)=norm_func(fun_k_plus_1_cell{1,i}(:),x_k_plus_1_cell{1,i}(:),spec.norm_spec(i));
        end


        mu_val_k_plus_1=mu_val_k;
        merit_val_vec_k_plus_1=merit_val_vec_k;
        
        if spec.line_search_spec==1

            %%% Implicitly assume common step size for all the variables %%%%

            obj_val_grad=compute_JVP_func(x_k_cell,d_k_cell,...
                obj_func,other_input_obj_func);%%%%

            d_norm=0;%L1 norm
            n_d=0;% number of elements in d_k_cell
            for i=1:size(x_k_cell,2)
                d_norm=d_norm+sum(abs(d_k_cell{i}(:)));
                n_d=n_d+size(d_k_cell{i}(:),1);
            end

            eta_k=0.1/((1+k)^2)*n_d;%%%%%

            kappa_val=0.5;
            c_val=0;%%%%
            numer=obj_val_grad+c_val*d_norm^2;
            denom=(1-kappa_val)*merit_val_vec_k(3)+eta_k;
            mu_val_k_plus_1=max(mu_val_k,numer/denom);
            
            %%%mu_val_k_plus_1=mu_val_k;%%%%%%%

            Delta_k=obj_val_grad-mu_val_k_plus_1*merit_val_vec_k(3);% Directional derivative of the merit function
            merit_val_vec_k_plus_1=compute_merit_val_func(x_k_plus_1_cell,mu_val_k_plus_1,...
                    obj_func,Phi_func,...
                    other_input_obj_func,other_input_Phi_func,...
                        spec.fixed_point_iter_spec,param_max,param_min);


            merit_val_k=merit_val_vec_k(2)+mu_val_k_plus_1*merit_val_vec_k(3);

            zeta_val=1e-4;
            if merit_val_vec_k_plus_1(1)>merit_val_k+...
                zeta_val*step_size*spectral_coef_k{1}(1)*Delta_k+...
                (1-zeta_val)*(1+mu_val_k_plus_1)*step_size*spectral_coef_k{1}(1)*eta_k
                    continue_backtracking_dummy=1;
            else
                continue_backtracking_dummy=0;    
            end

            rho=spec.rho;

            if continue_backtracking_dummy==1
                step_size=step_size.*rho;

            else
                break;
            end

        else% line search spec==0
            merit_val_vec_k_plus_1=compute_merit_val_func(x_k_plus_1_cell,mu_val_k_plus_1,...
                obj_func,Phi_func,...
                other_input_obj_func,other_input_Phi_func,...
                    spec.fixed_point_iter_spec,param_max,param_min);
        end%% line search spec==1 or 0

        if spec.line_search_spec==1 & iter_line_search==spec.ITER_MAX_LINE_SEARCH
            "Maximum number of line search reached."
        end
    end % end iter_line_search=1:ITER_MAX_LINE_SEARCH loop

    %merit_val_vec_k_plus_1
end
