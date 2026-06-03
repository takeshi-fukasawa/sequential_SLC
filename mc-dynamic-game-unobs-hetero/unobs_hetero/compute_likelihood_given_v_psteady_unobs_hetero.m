function [likelihood_vec,dobj_dparam,dobj_dv,dobj_dpsteady_1]=...
    compute_likelihood_given_v_psteady_unobs_hetero(...
    param,v,psteady_1,...
    middle_vars, numz_exo_vec,bdisc,...
        aval,mstate,dummy_aobs,dummy_xobs_obs,indobs_state_obs)

    %psteady_1=max(psteady_1,1e-16);%%%

pi_param=param(end);

n_param=size(param(:),1);
n_state=size(v,1);
n_state_unobs=numz_exo_vec(end);
n_state_obs=n_state/n_state_unobs;



[T,nobs,nplayer,nalt]=size(dummy_aobs);


    % Computing probabilities
    phat=reshape(v,[],nalt);
    phatmax = max(phat, [], 2);
    phat = phat - phatmax;
    numer=exp(phat);
    phat = numer ./ sum(numer, 2);%(n_state*nplayer)*nalt
    phat_reshaped=reshape(phat,n_state_obs,n_state_unobs,nplayer,nalt);


    phat_obs=sum(...
      phat_reshaped(indobs_state_obs(:),:,:,:).*...
      reshape(dummy_aobs,T*nobs,1,nplayer,nalt),...
      4);%(T*nobs)*n_state_unobs*nplayer*1

   prod_phat_obs=reshape(prod(phat_obs,3),T,nobs,n_state_unobs);

   [ptran_unobs,dptrans_unobs_dpi]=compute_ptrans_func(pi_param,n_state_unobs);
    
    if nargout==1
        [likelihood_vec]=  compute_Pr_x_a_func_unobs_hetero(prod_phat_obs,psteady_1,indobs_state_obs,ptran_unobs);
    
    else % Compute gradient
    
        [likelihood_vec,dl_dprod_phat_obs_array,dl_dpsteady_1,dl_dptrans_unobs]=...
            compute_Pr_x_a_func_unobs_hetero(prod_phat_obs,psteady_1,indobs_state_obs,ptran_unobs);
        
        dl_dpi=sum(dl_dptrans_unobs.*reshape(dptrans_unobs_dpi,1,n_state_unobs,n_state_unobs),[2,3]);%nobs*1
        dl_dparam=cat(2,zeros(nobs,n_param-1),dl_dpi);%nobs*nparam
        dobj_dparam_temp=-sum(dl_dparam./likelihood_vec,1)./nobs;
        dobj_dparam=dobj_dparam_temp(:);

        dobj_dpsteady_1=-sum(dl_dpsteady_1./likelihood_vec,1)./nobs;
        dobj_dpsteady_1=dobj_dpsteady_1(:);

        %%% dl_dprod_phat_obs_array: T*nobs*n_state_unobs
        
        dphat_a0_dv_a0=reshape(phat(:,1).*(1-phat(:,1)),n_state_obs,n_state_unobs,nplayer);
        dphat_a0_dv_a1=reshape(-phat(:,1).*phat(:,2),n_state_obs,n_state_unobs,nplayer);
        dphat_a1_dv_a0=reshape(-phat(:,1).*phat(:,2),n_state_obs,n_state_unobs,nplayer);
        dphat_a1_dv_a1=reshape(phat(:,1).*(1-phat(:,1)),n_state_obs,n_state_unobs,nplayer);
        
        dphat_a0_dv=cat(4,dphat_a0_dv_a0,dphat_a0_dv_a1);%n_state_obs*n_state_unobs*nplayer*nalt
        dphat_a1_dv=cat(4,dphat_a1_dv_a0,dphat_a1_dv_a1);%n_state_obs*n_state_unobs*nplayer*nalt
        dphat_dv=cat(5,dphat_a0_dv,dphat_a1_dv);%n_state_obs*n_state_unobs*nplayer*nalt(v)*nalt(phat)
        
        %%%%%%%%
        weight_temp=reshape(dl_dprod_phat_obs_array,T,nobs,n_state_unobs,1,1).*...
            reshape(prod_phat_obs,T,nobs,n_state_unobs,1,1)./...
            reshape(phat_obs,T,nobs,n_state_unobs,nplayer,1)./...
            reshape(likelihood_vec,1,nobs,1,1,1).*...
            reshape(dummy_aobs,T,nobs,1,nplayer,nalt);%T*nobs*n_state_unobs*nplayer*nalt; 
        weight=reshape(weight_temp,T*nobs,n_state_unobs,nplayer,nalt);
        
        weight_agg = zeros(n_state_obs, n_state_unobs, nplayer,nalt);
        for p = 1:nplayer
            for u = 1:n_state_unobs
                for m=1:nalt
                    weight_agg(:, u, p,m) = accumarray(indobs_state_obs, weight(:, u, p,m));% Aggregation concerning observed states
                end
            end
        end

        dobj_dv_temp=...
            -sum(reshape(weight_agg,n_state_obs,n_state_unobs,nplayer,1,nalt).*...
            reshape(dphat_dv,n_state_obs,n_state_unobs,nplayer,nalt,nalt),5);%n_state_obs*n_state_unobs*nplayer*nalt(v)
        dobj_dv=reshape(dobj_dv_temp,n_state*nplayer*nalt,1)./nobs;%(n_state*nplayer*nalt)*1
        dobj_dv=dobj_dv(:);
    end%nargout==2

end