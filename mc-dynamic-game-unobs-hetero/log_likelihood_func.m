function [obj_val,grad]=log_likelihood_func(b0,ydum,e0,e1,u0,u1,eobs0,eobs1,uobs0,uobs1)

global scale_ll


 x=[uobs0, uobs1];
 restx= [eobs0, eobs1];
 scale_ll=1;

   %%%% Negative values for minimization %%%%%
  myzero = 1e-16;
  nobs = size(ydum, 1);
  nalt = max(ydum);
  npar = size(x, 2) / nalt;


    % Computing probabilities
    phat = zeros(nobs, nalt);
    %%x_reshaped=reshape(x,nobs,npar,nalt);
    for j = 1:nalt
      phat(:,j) = x(:,npar*(j-1)+1:npar*j)*b0 + restx(:,j);
      %%phat(:,j) = x_reshaped(:,:,j)*b0 + restx(:,j);
    end
    phatmax = repmat(max(phat, [], 2), 1, nalt);
    phat = phat - phatmax;
    exp_phat=exp(phat);
    phat = exp_phat ./ repmat(sum(exp_phat, 2), 1, nalt);

    % Computing xmean
    llike = 0;%-nobs;
    for j = 1:nalt
      llike = llike ...
              + sum( (ydum==j) ...
                     .* log( (phat(:,j) > myzero).*phat(:,j) ...
                             + (phat(:,j) <= myzero).*myzero) );
    end

    llike=llike./nobs;%%%%

   obj_val=llike*(-1);%%%%%%%

   obj_val=obj_val./scale_ll;%%%%%%%%

    if nargout>=2
        % Computing gradient
        xysum = 0;
        for j = 1:nalt
            xysum = xysum + sum( repmat(ydum==j, 1, npar) .* x(:,npar*(j-1)+1:npar*j) );
        end
    
        sumpx = zeros(nobs, npar);
        for j=1:nalt
            xbuff = x(:,npar*(j-1)+1:npar*j);
            sumpx = sumpx + repmat(phat(:,j), 1, npar) .* xbuff;
        end
    
        d1llike = xysum - sum(sumpx);
    
        d1llike=d1llike./nobs;%%%%
        
        grad=d1llike*(-1);%%%%%%

        grad=grad/scale_ll;%%%%%

        grad=reshape(grad,size(b0));
    end% nargout>=2

end % function


