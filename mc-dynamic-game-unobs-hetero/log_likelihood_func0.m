function [obj_val,grad]=log_likelihood_func0(...
    b0,dummy_count,e0,e1,u0,u1)

 global scale_ll

 x=[u0, u1];
 restx= [e0, e1];

    cconvb = 1e-6;
  myzero = 1e-16;
  nobs = sum(dummy_count(:), 1);
  nalt = size(dummy_count,2);
  npar = size(x, 2) / nalt;
  err = 0;

 scale_ll=1;

   % x: (n_state*nPlayer)*(npar*nalt)

  xysum=sum(reshape(x,[],npar,nalt).*reshape(dummy_count,[],1,nalt),[1,3]);%1*npar

    % Computing probabilities

    phat=reshape(sum(reshape(x,[],npar,nalt).*reshape(b0,1,npar,1),2),[],nalt)+ restx;
    phatmax = max(phat, [], 2);
    phat = phat - phatmax;
    phat = exp(phat) ./ sum(exp(phat), 2);%[]*nalt


    llike=sum(log(phat).*dummy_count,[1,2]);

     llike=llike./nobs;

    obj_val=llike*(-1);%%%%%%%

   obj_val=obj_val./scale_ll;

    if nargout>=2
        % Computing gradient
          dummy_count_sum=sum(dummy_count,2);

          mean_x_temp=sum(reshape(phat,[],1,nalt).*reshape(x,[],npar,nalt),3);%[]*npar
    
          mean_x=sum(reshape(mean_x_temp,[],npar).*...
            reshape(dummy_count_sum,[],1),1);%1*npar
    
        d1llike=xysum-mean_x;
        d1llike=d1llike./nobs;
        grad=d1llike*(-1);%%%%%%
        grad=reshape(grad,size(b0));
    end

end % function


