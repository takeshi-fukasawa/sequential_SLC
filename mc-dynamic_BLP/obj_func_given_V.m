function [obj,grad]=obj_func_given_V(param,V,...
weight,v,X_jt,X,Z,PX,W,...
        S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec)

    mu_ijt=mu_func2(param,v,X_jt);

    [J,ns,G,T]=size(mu_ijt);

    n_dim_V=size(V,4);

    if n_dim_V>T %% IVS spec
       V=V(:,:,:,1:T);%1*ns*1*T; V at obs data pts
    end

    s_ijt_ccp_up_to_scale=exp(reshape(mu_ijt,J,ns,G,T)-...
        reshape(V,1,ns,1,T));%J*ns*G*T


   delta=NaN(J,1,G,T);
   Pr0=ones(1,ns,1,T);
    for t=1:T
           Pr0_t=Pr0(:,:,:,t,:);
           s_jt_up_to_scale_t=...
               sum(reshape(weight,1,ns,1,1).*Pr0_t.*...
               s_ijt_ccp_up_to_scale(:,:,:,t),2);%J*1*G*1
    
          exp_delta_t=S_jt_data(:,:,:,t)./s_jt_up_to_scale_t;%J*1*G*1
          delta(:,:,:,t)=log(exp_delta_t);%J*1*G*1
      
          s_ijt_ccp_t=s_ijt_ccp_up_to_scale(:,:,:,t,:).*...
              exp_delta_t;%J*ns*G*1

          s_i0t_ccp_t=1-sum(s_ijt_ccp_t,[1,3]);%1*ns*1*1

          %%% Avoid complex numbers
          small_val=1e-8;
          s_i0t_ccp_t=max(s_i0t_ccp_t,small_val);%1*ns*1*1

          if t<=T-1
             Pr0(:,:,:,t+1)=Pr0_t.*s_i0t_ccp_t;%1*ns*1*1
          end
    end % for loop


    obj=obj_func_given_delta(delta,PX,X,Z,W);

    grad=[];
end

