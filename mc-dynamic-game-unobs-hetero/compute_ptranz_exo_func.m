function [ptranz_exo,grad]=compute_ptranz_exo_func(pi_param,numz_exo_vec)

    if size(numz_exo_vec(:),1)==1
        if nargout==1
            ptranz_exo=compute_ptrans_func(pi_param,numz_exo_vec);
        else
            [ptranz_exo,grad]=compute_ptrans_func(pi_param,numz_exo_vec);
        end  
    elseif size(numz_exo_vec(:),1)==2

        n_state_obs=numz_exo_vec(1);
        n_state_unobs=numz_exo_vec(2);


        if nargout==1
            ptranz_exo_unobs=compute_ptrans_func(pi_param,n_state_unobs);
        else
            [ptranz_exo_unobs,grad_unobs]=compute_ptrans_func(pi_param,n_state_unobs);
        end        

        pi_param_obs=0.8;
        ptranz_exo_obs=compute_ptrans_func(pi_param_obs,n_state_obs);

        ptranz_exo=kron(ptranz_exo_unobs, ones(n_state_obs,n_state_obs)).*...
            repmat(ptranz_exo_obs,n_state_unobs,n_state_unobs);
        %(n_state_obs*n_state_unob)*(n_state_obs*n_state_unob)

        if nargout==2
            grad=kron(grad_unobs, ones(n_state_obs,n_state_obs)).*...
                repmat(ptranz_exo_obs,n_state_unobs,n_state_unobs);
            %(n_state_obs*n_state_unob)*(n_state_obs*n_state_unob)
        end
    end%size(numz_exo_vec(:),1)==2 or 1

end
