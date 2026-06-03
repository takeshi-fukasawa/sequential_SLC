function mu_ijt=mu_func2(param,v,X_jt)

    global positive_param_nfxp_spec
    if positive_param_nfxp_spec==1
        param=abs(param);% Make the std of random coeffcients positive
    end
    
    sigma_vec=[0,param(1),param(2),0,param(3)];
   [J,~,G,T,K]=size(X_jt);
   ns=size(v,2);
   mu_ijt=sum(reshape(sigma_vec,1,1,1,1,K).*...
       reshape(v,1,ns,1,1,K).*...
       reshape(X_jt,J,1,G,T,K),5);%J*ns*G*T
end
