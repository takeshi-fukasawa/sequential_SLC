function [output,other_vars]=V_param_update_func(V_initial,param,...
        v,X_jt,weight,S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec)


[J,~,G,T,~]=size(X_jt);
ns=size(weight,2);


mu_ijt=mu_func2(param,v,X_jt);

    [output,other_vars]=V_update_func(...
        V_initial,weight,mu_ijt,...
        S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec);

    if 1==0 % The following works well in rational expectation (perfect foresight) model
        %%% Essentially nonstationary environment %%%

        %mat=triu(ones(T))';%T*T
        mat=tril(ones(T));%T*T
        
        V_updated=output{1};
        resid=V_updated(:,:,:,1:T)-V_initial(:,:,:,1:T);%1*ns*1*T
        resid_new=sum(reshape(resid,1,ns,1,T,1).*...
            reshape(mat,1,1,1,T,T),4);%1*ns*1*1*T
    
        V_updated(:,:,:,1:T)=V_initial(:,:,:,1:T)+...
            reshape(resid_new,1,ns,1,T);
        output{1}=V_updated;
    end
    
end