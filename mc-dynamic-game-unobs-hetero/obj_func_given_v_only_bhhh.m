function ll_i=obj_func_given_v_bhhh(param,v,middle_vars,... 
 ptranz_exo,bdisc,...
        aval,mstate,dummy_count,ptranz_exo_kron)


  nobs = sum(dummy_count(:), 1);
  nalt = size(dummy_count,2);

    % Computing probabilities

    phat=reshape(v,[],nalt);
    phatmax = max(phat, [], 2);
    phat = phat - phatmax;
    phat = exp(phat) ./ sum(exp(phat), 2);%[]*nalt

    ll_i=sum(log(phat).*dummy_count,2);%n_obs*1

end
