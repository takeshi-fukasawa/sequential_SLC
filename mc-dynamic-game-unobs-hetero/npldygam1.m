% ----------------------------------------------------------------------------
% NPLDYGAM       Procedure that estimates the structural parameters
%                of dynamic game of firms' entry/exit
%                using a Nested Pseudo-Likelihood (NPL) algorithm
%
% Based on the NPL Gauss code by Victor Aguirregabiria.
% Converted from Gauss to Matlab by Jason Blevins.
% Modified to include CCPs and (un-normalized) choice-specific value
% functions from 1-step NPL along with return values.
%
% ----------------------------------------------------------------------------
%
% FORMAT:
% [bnpl,varb,pchoice,v1npl,vnpl,iter,err,times] = npldygam(aobs,zexo_obs,aobs_1,z_exoval,ptranz_exo,pchoice,bdisc,miniter,maxiter)
%
% INPUTS:
%     aobs    - (nobs x nplayer) matrix with observations of
%               firms' activity decisions (1=active ; 0=no active)
%     zexo_obs    - (nobs x 1) vector with observations of market
%               exogenous characteristics.
%     aobs_1  - (nobs x nplayer) matrix with observations of
%               firms' initial states (1=incumbent; 0=potential entrant)
%     z_exoval    - (numz_exo x 1) vector with values of market characteristics
%     ptranz_exo  - (numz_exo x numz_exo) matrix of transition probabilities
%               of market characteristics.
%     pchoice - ((numz*2^nplayer) x nplayer) matrix of players'
%               choice probabilities used to initialize the procedure.
%     bdisc   - Discount factor
%     miniter - Minimum number of NPL iterations
%     maxiter - Maximum number of NPL iterations
%
% OUTPUTS:
%     bnpl    - (kparam x niter) matrix with NPL estimates as rows:
%               (theta_fc_1; theta_fc_2; ... ; theta_fc_n; theta_rs; theta_rn; theta_ec).
%     varb    - (kparam x kparam) variance matrix from clogit estimation.
%     pchoice - ((numz*2^nplayer) x nplayer) matrix of players'
%               choice probabilities consistent with estimates.
%     v1npl   - (nplayer x 2 x (numz*2^nplayer)) matrix of players'
%               choice-specific value functions consistent with estimates
%               from the first iteration of NPL.
%     vnpl    - (nplayer x 2 x (numz*2^nplayer)) matrix of players'
%               choice-specific value functions consistent with the
%               converged NPL estimates.
%     iter    - Number of NPL iterations completed.
%     err     - Error flag indicator
%     times   - Vector of computational time (sec.) per iteration
%
% ----------------------------------------------------------------------------

function [ bnpl, varb, pchoice, v1npl, vnpl, llike, iter, err, times ] = ...
    npldygam1(aobs, zexo_obs, aobs_1, z_exoval, ptranz_exo, pchoice, bdisc, miniter, maxiter)

  run preliminary_construct_variables.m

bcnpl = zeros(kparam,1);
bnpl = [];
  v1npl = zeros(n_state, nplayer, naction);
  vnpl = zeros(n_state, nplayer, naction);

  indobs=indobs_orig;

  % -------------
  % NPL algorithm
  % -------------
  aobs = 1 + reshape(aobs, [nobs*nplayer,1]);
  v = zeros(n_state,nplayer,naction);

 
  % Iterate until convergence criterion met
  criter = 1000;
  iter = 1;
  while ((criter > critconv) && (iter <= maxiter))
    %fprintf('\n');
    fprintf('-----------------------------------------------------\n');
    fprintf('NPL ESTIMATOR: K = %d\n', iter);
    fprintf('-----------------------------------------------------\n');
    fprintf('\n');
    tic;


    out_e_u=...
        compute_e_u_from_pchoice_func(pchoice, ptranz_exo,  bdisc,...
        aval,mstate,ptranz_exo_kron);

    e0=out_e_u{1};
    e1=out_e_u{2};
    u0=out_e_u{3};
    u1=out_e_u{4};

    %%% eobs0,eobs1,uobs0,uobs1 %%%%
    % ------------------------------------------
    %  (d) Pseudo Maximum Likelihood Estimation
    % ------------------------------------------
    %tic
    %%[ thetaest, varb, llike, cl_iter, err ] = clogit(aobs, [uobs0, uobs1], [eobs0, eobs1], 100, bcnpl);
    [ thetaest, varb, llike, cl_iter, err ] = clogit0(dummy_count,[u0,u1],[e0,e1],...
        100, bcnpl);

    %toc

    %if err > 0
    %    return;
    %end



   %err=0;%%%%%%
    %%%%%%varb


    % ----------------------------
    %  (e) Updating probabilities
    % ----------------------------
    for i = 1:nplayer
      v(:,i,2) = u1((i-1)*n_state+1:i*n_state,:)*thetaest + e1((i-1)*n_state+1:i*n_state,:);
      v(:,i,1) = u0((i-1)*n_state+1:i*n_state,:)*thetaest + e0((i-1)*n_state+1:i*n_state,:);
      buff = v(:,i,2) - v(:,i,1);
      pchoice_est(:,i) = exp(buff)./(1+exp(buff));
    end

    % Check for convergence after minimum number of iterations
    if (iter > miniter)
        criter1 = max(abs(thetaest - bcnpl));
        criter2 = max(max(abs(pchoice_est - pchoice)));
        criter = max(criter1, criter2);
    end

    % Save values for next iteration
    bcnpl = thetaest;
    vnpl = v;
    pchoice = pchoice_est;

    % Collect estimates
    bnpl = [ bnpl; thetaest' ];

    % Collect times
    times = [ times; toc ];

    % Separately save value function from first iteration
    if iter == 1
        v1npl = v;
    end

    disp('theta = ')
    disp(bcnpl')
    disp('llike = ')
    disp(llike)

    % Proceed to the next iteration
    iter = iter + 1;
  end

  % Undo final iteration increment
  iter = iter - 1;
end
