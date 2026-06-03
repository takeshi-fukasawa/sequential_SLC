function [Pr_x_a_vec,dl_dprob_array,dl_dpsteady_1,dl_dptran_unobs]=...
    compute_Pr_x_a_func_unobs_hetero(prod_phat_obs,psteady_1,indobs_state_obs,ptran_unobs)

[T,nobs,n_state_unobs]=size(prod_phat_obs);
n_state=size(psteady_1(:),1)+1;
n_state_obs=n_state/n_state_unobs;

if nargout==1
  [Pr_s1_given_x1]=compute_Pr_s1_given_x1_func(psteady_1,n_state_obs,n_state_unobs);
else % Compute Jacobians
  %[Pr_s1_given_x1,Jac_dPr_s1_given_x1_dpsteady_1]=...
  %    compute_Pr_s1_given_x1_func(psteady_1,n_state_obs,n_state_unobs);
  [Pr_s1_given_x1,Jac_dPr_s1_given_x1_dpsteady_1]=...
      compute_Pr_s1_given_x1_func2(psteady_1,n_state_obs,n_state_unobs);

end

indobs_state_obs_reshaped=reshape(indobs_state_obs,T,nobs);
Pr_s1_given_x1_obs=Pr_s1_given_x1(reshape(indobs_state_obs_reshaped(1,:),[],1),:);%n_obs*n_state_unobs

if nargout==1
  h_t=Pr_s1_given_x1_obs;
else
  h_t_array=NaN(T,nobs,n_state_unobs);
  h_t_array(1,:,:)=Pr_s1_given_x1_obs;%nobs*n_state_unobs%%%
end

  %% Sequentially compute h_t
for t=1:T

   prod_phat_obs_t=reshape(prod_phat_obs(t,:,:),nobs,n_state_unobs);

   if nargout>=2
     h_t=reshape(h_t_array(t,:,:),nobs,n_state_unobs);
   end

    if t==T
       Pr_x_a_vec=sum(prod_phat_obs_t.*h_t,2);%nobs*1
   end

  prob_at_xt1_given_xt_st=prod_phat_obs_t;%nobs*n_unobs

  h_t1=sum(reshape(ptran_unobs,1,n_state_unobs,n_state_unobs).*...
      reshape(prob_at_xt1_given_xt_st.*h_t,nobs,n_state_unobs,1),2);%nobs*1*n_state_unobs

  if nargout==1
    h_t=reshape(h_t1,nobs,n_state_unobs);
  else
    h_t_array(t+1,:,:)=reshape(h_t1,nobs,n_state_unobs);%%%
  end

end % for loop wrt t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Jacobian
if nargout>=2
  dht1_dht_array=reshape(ptran_unobs,1,1,n_state_unobs,n_state_unobs).*...
  reshape(prod_phat_obs,T,nobs,n_state_unobs,1);%T*nobs*n_state_unobs*n_state_unobs

  dl_dh_array=NaN(T,nobs,n_state_unobs);
  dl_dh_array(T,:,:)=reshape(prod_phat_obs(T,:,:),nobs,n_state_unobs);

  dl_dprob_array=NaN(T,nobs,n_state_unobs); % Derivative of likelihood wrt  prob_at_xt1_given_xt_st(t=1,...,T)
  dl_dprob_array(T,:,:)=reshape(h_t_array(T,:,:),nobs,n_state_unobs);

  dl_dptran_unobs_array=zeros(T,nobs,n_state_unobs,n_state_unobs); % Derivative of likelihood wrt  ptran_unobs(s_{t+1}|s_{t}; t=1,...,T-1)

  for t=T-1:-1:1
    dl_dh_array(t,:,:)=sum(...
    reshape(dl_dh_array(t+1,:,:),nobs,1,n_state_unobs).*...
    reshape(dht1_dht_array(t,:,:),nobs,n_state_unobs,n_state_unobs),3);%nobs*n_state_unobs(at time t)*1

    dl_dprob_array(t,:,:)=sum(...
    reshape(dl_dh_array(t+1,:,:),nobs,1,n_state_unobs).*...
    reshape(ptran_unobs,1,n_state_unobs,n_state_unobs).*...
    reshape(h_t_array(t,:,:),nobs,n_state_unobs,1),3);%nobs*n_state_unobs(at time t)
    
    dl_dptran_unobs_array(t,:,:,:)=...
     reshape(dl_dh_array(t+1,:,:),nobs,1,n_state_unobs).*...
     reshape(prod_phat_obs(t,:,:),nobs,n_state_unobs,1).*...
     reshape(h_t_array(t,:,:),nobs,n_state_unobs,1);%nobs*n_state_unobs(at time t)*n_state_unobs(at time t+1)

  end%fot t=T-1:-1:1

  dl_dPrs1_given_x1=dl_dh_array(1,:,:);

  Jac_dPr_s1_given_x1_dpsteady_1=reshape(Jac_dPr_s1_given_x1_dpsteady_1,n_state_obs,n_state_unobs,n_state-1);
  Jac_dPr_s1_given_x1_dpsteady_1_obs=...
       Jac_dPr_s1_given_x1_dpsteady_1(reshape(indobs_state_obs_reshaped(1,:),[],1),:,:);%n_obs*n_state_unobs*(n_state-1)
   
   dl_dpsteady_1_temp=reshape(dl_dPrs1_given_x1,nobs,n_state_unobs,1).*...
        reshape(Jac_dPr_s1_given_x1_dpsteady_1_obs,nobs,n_state_unobs,n_state-1);%nobs*n_state_unobs*(n_state-1)
   dl_dpsteady_1=reshape(sum(dl_dpsteady_1_temp,2),nobs,1,n_state-1);

   dl_dptran_unobs=reshape(sum(dl_dptran_unobs_array,1),nobs,n_state_unobs,n_state_unobs);
end

end% function
