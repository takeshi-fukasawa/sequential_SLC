%%%%%%%%%
function [v0SP,prob0SP,best_logit]=estimate_initial_CCP_func(...
    aobs,aobs_1,zexo_obs,mstate,ignore_zexo2_spec,...
    nobs,nplayer,nstate,numz_exo_vec)

naction=2;

  % Construct dependent (aobsSP) and explanatory variables (xobsSP)
  aobsSP = reshape(aobs', nobs*nplayer, 1);
  alphai = kron(ones(nobs,1), eye(nplayer));

    if ignore_zexo2_spec==0
      xobsSP = kron(zexo_obs, ones(nplayer,1));
      mstateSP = [ ones(size(mstate,1), nplayer), mstate, ...
                sum(mstate(:,2:nplayer+1)')' ];%nstate*n_vars

    elseif ignore_zexo2_spec==1
       xobsSP = kron(zexo_obs(:,1), ones(nplayer,1));
       mstateSP = [ ones(size(mstate,1), nplayer), mstate(:,1:end-1), ...
                sum(mstate(:,2:nplayer+1)')' ];%nstate*n_vars
       numz_exo_vec=numz_exo_vec(1);
    end

  nfirms_1 = kron(sum(aobs_1, 2), ones(nplayer,1));
  aobs_1SP = reshape(aobs_1', nobs*nplayer, 1);
  xobsSP = [ alphai, xobsSP, aobs_1SP, nfirms_1 ];
  % Logit estimation
  [ best_logit, varest ] = milogit(aobsSP, xobsSP);
  % Construct probabilities


  if size(numz_exo_vec,2)==1
    best = [ diag(best_logit(1:nplayer)) ; ...
            ones(1,nplayer) * best_logit(nplayer+1); ...%zexo
            eye(nplayer) * best_logit(nplayer+2); ...
            ones(1,nplayer) * best_logit(nplayer+3) ];
  else % Multiple exogenous state variables
    best = [ diag(best_logit(1:nplayer)) ; ...
            ones(1,nplayer) * best_logit(nplayer+1); ...% first zexo vars
            eye(nplayer) * best_logit(nplayer+3); ...
            ones(1,nplayer) * best_logit(nplayer+2); ...% second zexo vars
            ones(1,nplayer) * best_logit(nplayer+4) ];
  end


  %mstateSP:nstate*n_vars
  % best:n_vars*nplayer
  v0SP = zeros(nstate, nplayer, naction);
  v0SP(:,:,2) = mstateSP*best;
  prob0SP = 1 ./ (1+exp(-mstateSP*best));

  if ignore_zexo2_spec==1
      best_logit=[best_logit(1:nplayer+1);0;best_logit(nplayer+2:end)];% Let the coef of unobs state variable be zero
  end

end
