% mc_sequential.m 
% This is based on a Matlab translation of the Gauss program of
% Aguirregabiria and Mira (2007) originally named
% am_econometrica_2007_montecarlo.prg, with the EPL estimator added.
% Comments from the original are preserved below.
%
% Adam Dearing and Jason Blevins
% Columbus, Ohio
% June 30, 2020

% Modified by Takeshi Fukasawa

% ****************************************************************************
%
%  am_econometrica_2007_montecarlo.prg
%
%      GAUSS program that implements the Monte Carlo experiments reported in
%      Section 4 of the paper: Aguirregabiria, V., and P. Mira (2007):
%      'Sequential Estimation of Dynamic Discrete Games,' Econometrica,
%      Vol. 75, No. 1, pages 1-53.
%
%  by Victor Aguirregabiria
%
%  - Last version: May 2005
%  - Comments and explanations added on August 2011
%
% ****************************************************************************


% ------------------------------------------------------------------------------
%                              MODEL
% ------------------------------------------------------------------------------
%
%  MAIN FEATURES OF THE MODEL
%
%  - Dynamic game of firm entry-exit in a market
%
%  -   i = Firm index, i belongs to (1,2,3,4,5)
%      t = Time index
%
%  -   Decision variable:
%      a[it] = Indicator of the event 'firm i operates in the market at period t'
%
%          a[it] = 0 ----> No active in market
%          a[it] = 1 ----> Active in market
%
%  -   The game is dynamic because there is a sunk cost of entry
%
%  -   The state variables of this game are the indicators of whether
%      the firms were active or not in the market at previous period,
%      and market size. These are payoff relevant state variables because
%      they determine whether a firm has to pay an entry cost to operate
%      in the market.
%
%  -   State variables
%      x[t] = ( s[t], a[1,t-1], a[2,t-1], a[3,t-1], a[4,t-1], a[5,t-1] )
%
%      e0[it] and e1[it] = Firm i's private information.
%
%  - Profit function
%
%      If a[it] = 0 (firm i is not active at period t)
%
%          Profit(0) = e0[it]
%
%      If a[it] = 1 (firm i is not active at period t)
%
%          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1])
%                    + theta_rs * s[t]
%                    - theta_rn * ln(N[-it]) + 1) + e1[it]
%  where:
%          theta_fc_i, theta_ec, theta_rs, theta_rn are parameters
%
%          theta_fc_i  =   Firm i's fixed effect
%          theta_ec    =   Entry cost
%          N[-it]      =   Number of firms, other than i, operating in the market at period t.
%
%  eps_i(0) and eps_i(1) are private information shocks that are
%  i.i.d. over time, markets, and players, and independent of each other,
%  with Extreme Value Type 1 distribution.
%
%
% ------------------------------------------------------------------------------
%              MONTE CARLO EXPERIMENTS
% ------------------------------------------------------------------------------

% -----------------------------------------------------------
%  SELECTION OF THE MONTE CARLO EXPERIMENT TO IMPLEMENT
% -----------------------------------------------------------
global compute_jacobian_spec krylov_spec

if ~exist('selexper')
    selexper = 1; % Value from 1 to 6 that represents the index of the
                  % Monte Carlo experiment to implement. A run of this
                  % program implements one experiment.
end

if ~exist('nobs')
    nobs = 1600;  % Number of markets (observations)
end

% -----------------------------------------------
%  VALUES OF PARAMETERS AND OTHER CONSTANTS
% -----------------------------------------------
naction = 2;    % Number of actions
numexp = 3;     %  Total number of Monte Carlo experiments
nstart = 1;     % Number of starting values
T=10; % Observations in each market
miniter = 3;    %  Minimum number of NPL or EPL iterations
maxiter = 100; %  Maximum number of NPL or EPL iterations
maxiter_npl=3;

theta_fc = zeros(numexp, nplayer);
%theta_fc(:, 1) = -1.9 * ones(numexp, 1); % Vector with values of theta_fc_1 for each experiment
%theta_fc(:, 2) = -1.8 * ones(numexp, 1); % Vector with values of theta_fc_2 for each experiment
%theta_fc(:, 3) = -1.7 * ones(numexp, 1); % Vector with values of theta_fc_3 for each experiment
%theta_fc(:, 4) = -1.6 * ones(numexp, 1); % Vector with values of theta_fc_4 for each experiment
%theta_fc(:, 5) = -1.5 * ones(numexp, 1); % Vector with values of theta_fc_1 for each experiment

theta_fc(:, :) = repmat(-2.0+0.1*[1:nplayer],numexp,1); % Vector with values of theta_fc_1 for each experiment


theta_rs = 1.0 * ones(numexp, 1); % Vector with values of theta_rs for each experiment

theta_rs2 = 1.0 * ones(numexp, 1); % Vector with values of theta_rs2 for each experiment

disfact = 0.95 * ones(numexp, 1); % Vector with values of discount factor for each experiment
sigmaeps = 1 * ones(numexp, 1);   % Vector with values of std. dev. epsilon for each experiment

theta_rn = [ 1.0; 2.5; 4.0 ]; % Vector with values of theta_rn for each experiment
theta_ec = [ 1.0; 1.0; 1.0 ]; % Vector with values of theta_ec for each experiment

% Points of support and transition probability of state variable s[t], market size
z_exo_val_orig = [ 1:1:5 ]';  % Support of market size
%z_exo_val_orig = [ 1:0.5:5]';  % Support of market size
z_exo2_val_orig= [ 1:1:3 ]';  % Support of market size (second variable)

z_exo_val=z_exo_val_orig;

%%%%
if multiple_zexo_spec==0
   numz_exo_vec=size(z_exo_val,1);
else
   numz_exo_vec=[size(z_exo_val_orig,1),size(z_exo2_val_orig,1)];
   z_exo_val=[repmat(z_exo_val_orig,size(z_exo2_val_orig,1),1),...
       kron(z_exo2_val_orig,ones(size(z_exo_val_orig,1),1))];
end

%%%%

numz_exo = size(z_exo_val, 1);  % Number of possible market sizes
nstate = numz_exo * (2^nplayer);  % Number of points in the state space
%ptranz_exo = [ 0.8, 0.2, 0.0, 0.0, 0.0 ;
%           0.2, 0.6, 0.2, 0.0, 0.0 ;
%           0.0, 0.2, 0.6, 0.2, 0.0 ;
%           0.0, 0.0, 0.2, 0.6, 0.2 ;
%           0.0, 0.0, 0.0, 0.2, 0.8 ]; % numz_exo=5 case

pi_param_true=0.8;
ptranz_exo=compute_ptranz_exo_func(pi_param_true,numz_exo_vec);


% Selecting the parameters for the experiment
theta_fc = theta_fc(selexper, :)';
theta_rs = theta_rs(selexper);
theta_rs2 = theta_rs2(selexper);
disfact = disfact(selexper);
sigmaeps = sigmaeps(selexper);
theta_ec = theta_ec(selexper);
theta_rn = theta_rn(selexper);

% Vector with true values of parameters
if size(numz_exo_vec(:),1)==1
   trueparam = [ theta_fc; theta_rs; theta_rn; theta_ec; disfact; sigmaeps ];
else
   trueparam = [ theta_fc; theta_rs;theta_rs2; theta_rn; theta_ec; disfact; sigmaeps ];
end



% Number of parameters to estimate
kparam = size(trueparam, 1) - 2;

% Structure for storing true parameters and settings
param.theta_fc = theta_fc;
param.theta_rs = theta_rs;
param.theta_rs2 = theta_rs2;
param.theta_rn = theta_rn;
param.theta_ec = theta_ec;
param.disfact = disfact;
param.sigmaeps = sigmaeps;
param.z_exo_val = z_exo_val;
param.ptranz_exo = ptranz_exo;
param.verbose = 1;

% Seed for (pseudo) random number generation
rand('seed', 20150403);

if 1==0% Misspecified model %%%
    disfact_est=0.99;%disfact;%%%
    pi_param_est=0.9;%pi_param_true*0.9;%%%%
    ptranz_exo_est=compute_ptranz_exo_func(pi_param_est,numz_exo_vec);%%%%%
else
    disfact_est=disfact;
    pi_param_est=pi_param_true;
    ptranz_exo_est=ptranz_exo;
end



numa = 2^nplayer;
aval = zeros(numa, nplayer);
for i = 1:nplayer
  aval(:,i) = kron(kron(ones(2^(i-1),1), [ 0; 1 ]), ones(2^(nplayer - i), 1));
end

mstate=compute_mstate_func(param.z_exo_val,nstate,nplayer,numa,numz_exo,aval);
ptranz_exo_kron=(kron(ptranz_exo, ones(numa,numa)));%%%

% -----------------------------------------------------------------------
%  COMPUTING A MARKOV PERFECT EQUILIBRIUM OF THIS DYNAMIC GAME
% -----------------------------------------------------------------------
eqfile = sprintf('mc_epl_eq%d.mat', selexper);
%if isfile(eqfile)
%    disp('Loading pre-calculated equilibrium from file.')
%    load(eqfile);
%else
    disp('Solving for an equilibrium of the game.')
    v0 = zeros(nstate*nplayer*naction,1);
    

    theta=trueparam(1:end-2);
    spec=[];
    spec.TOL=1e-8;
    %spec.update_spec=0;
    spec.Anderson_acceleration=1;
    [out,other_vars,iter_info_v_solve]=spectral_func(@v_update_func,spec,{zeros(nstate,nplayer,naction)},...
       theta,[],ptranz_exo,disfact,aval,mstate,ptranz_exo_kron);
    vequil=out{1};

    %eqcond = @(v) Gfunc(z_exo_val, ptranz_exo, trueparam(1:kparam), v, disfact, nplayer, naction);
    %vequil = fsolve(eqcond, v0, optimoptions('fsolve', 'Display', 'iter'));


    vmat = reshape(vequil, [nstate, nplayer, naction]);

    pequil = exp(vmat(:,:,2))./(exp(vmat(:,:,1))+exp(vmat(:,:,2)));

    %if save_spec==1
    %    save(eqfile, 'vequil', 'vmat', 'pequil');
    %end

%end

% -----------------------------------------------------------------------
%  COMPUTING THE STEADY STATE DISTRIBUTION OF THIS DYNAMIC GAME
% -----------------------------------------------------------------------
ptrana = ones(nstate, numa);
for i = 1:nplayer;
  mi = aval(:,i)';
  ppi = pequil(:,i);
  ppi1 = repmat(ppi, 1, numa);
  ppi0 = 1 - ppi1;
  mi1 = repmat(mi, nstate, 1);
  mi0 = 1 - mi1;
  ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
end
ptrana = (kron(param.ptranz_exo, ones(numa, numa))) .* (kron(ones(1,numz_exo), ptrana));
critconv = (1e-3)*(1/nstate);
criter = 1000;
psteady0 = (1/nstate) * ones(nstate, 1);
while (criter > critconv)
  psteady1 = ptrana' * psteady0;
  criter = max(abs(psteady1 - psteady0), [], 1);
  psteady0 = psteady1;
end
psteady = psteady0;

% --------------------------------------------------------------------------
%  SIMULATING DATA (50,000 OBSERVATIONS) FROM THE EQUILIBRIUM
%  TO OBTAIN DESCRIPTIVE STATISTICS ON THE DYNAMICS OF MARKET STRUCTURE
% --------------------------------------------------------------------------
nobsfordes = 50000;
[ a_mpe_test, a_mpe_test_1 ] = simdygam(nobsfordes, pequil, psteady, mstate);


mpestat(a_mpe_test, a_mpe_test_1);

% --------------------------------------------------------------------------
%  MONTE CARLO EXPERIMENT
% --------------------------------------------------------------------------

fprintf('\n');
fprintf('*****************************************************************************************\n');
fprintf('       MONTE CARLO EXPERIMENT #%d\n', selexper);
fprintf('*****************************************************************************************\n');
fprintf('\n');
param.verbose = 1;

% Cell arrays of starting values for each iteration
start_name = {}; % Name
start_name_short = {}; % Short name for starting value
start_p = {}; % CCPs
start_theta = zeros(kparam, nstart); % Parameters
start_v = {}; % Value function
start_1npl = zeros(nstart, 1); % Whether or not to use 1-NPL values to initialize EPL
start_err = zeros(nstart, 1); % Failure
start_iter = zeros(nstart, 1); % Iteration count
start_fail = zeros(nstart, 1); % Failure
start_ll = zeros(nstart, 1);
start_theta1 = zeros(kparam, nstart); % 1-step estimates for each starting value
start_thetac = zeros(kparam, nstart); % Converged estimates for each starting value
start_v1npl = {}; % 1-NPL value function for each starting value
theta_npl = {}; % Sequences of estimates for each starting value
theta_epl = {}; % Sequences of estimates for each starting value
theta_seq = {}; % Sequences of estimates for each starting value
theta_nfxp = {}; % Sequences of estimates for each starting value


start_times = {}; % Sequence of computational times per estimate

% Matrices of estimates for each iteration
bmat_1npl = zeros(nrepli, kparam); % 1-NPL estimates
bmat_2npl = zeros(nrepli, kparam); % 2-NPL estimates
bmat_3npl = zeros(nrepli, kparam); % 3-NPL estimates
bmat_cnpl = zeros(nrepli, kparam); % Converged NPL estimates
bmat_1epl = zeros(nrepli, kparam); % 1-EPL estimates
bmat_2epl = zeros(nrepli, kparam); % 2-EPL estimates
bmat_3epl = zeros(nrepli, kparam); % 3-EPL estimates
bmat_cepl = zeros(nrepli, kparam); % Converged EPL estimates
btrue_npl = zeros(nrepli, kparam); % 1-NPL estimates using true CCPs
btrue_epl = zeros(nrepli, kparam); % 1-EPL estimates using true value function and parameters

% Iteration counts for each method
iter_cnpl = zeros(nrepli, 1);
iter_cepl = zeros(nrepli, 1);

% Total computational times for each method
time_1npl = zeros(nrepli, 1);
time_1epl = zeros(nrepli, 1);
time_2npl = zeros(nrepli, 1);
time_2epl = zeros(nrepli, 1);
time_3npl = zeros(nrepli, 1);
time_3epl = zeros(nrepli, 1);
time_cnpl = zeros(nrepli, 1);
time_cepl = zeros(nrepli, 1);

% Per iteration times for each method
time_2npl_iter = zeros(nrepli, 1);
time_2epl_iter = zeros(nrepli, 1);
time_3npl_iter = zeros(nrepli, 1);
time_3epl_iter = zeros(nrepli, 1);
time_cnpl_iter = zeros(nrepli, 1);
time_cepl_iter = zeros(nrepli, 1);

% Failure counts for each method
fail_cnpl = zeros(nrepli, 1);
fail_cepl = zeros(nrepli, 1);

% Note: When there are multicollinearity problems in a Monte Carlo
% sample we ignore that sample and take a new one. We want to check
% for the number of times we have to make these redraws.
redraws = 0;

for draw = 1:nrepli
  fprintf('=========================================================================================\n');
  fprintf('     Replication = %d\n', draw);
  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.0)   Simulations of x''s and a''s\n');
  flag = 1;
  while (flag==1)
    [ aobs, aobs_1] = simdygam(nobs, pequil, psteady, mstate);
    [ aobs, aobs_1, zexo_obs] = simdygam_forward(nobs/T, pequil, psteady, mstate,ptrana,T);

    check = sum(sum([ aobs, aobs_1 ]) == zeros(1, 2*nplayer));
    check = check + sum(sum([ aobs, aobs_1 ]) == (nobs .* ones(1, 2*nplayer)));
    if (check > 0)
        flag = 1;
    elseif (check == 0)
        flag = 0;
    end
    redraws = redraws + flag; % Counts the number re-drawings
  end

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.1)   Estimation of initial CCPs (Semi-Parametric: Logit)\n');

%%%%%
ignore_zexo2_spec=0;%%% Even when zexo2 is unobserved, treat the variable as if observed.
[v0SP,prob0SP,best_logit]=estimate_initial_CCP_func(aobs,aobs_1,zexo_obs,mstate,ignore_zexo2_spec,...
    nobs,nplayer,nstate,numz_exo_vec);

  start_name{1} = 'Semi-Parametric: Logit';
  start_name_short{1} = 'SP';
  start_p{1} = prob0SP;
  start_1npl(1) = 1;
  start_theta(:,1) = best_logit;
  start_v{1} = v0SP;


  %--------------------------------------------------
  % NPL
  %--------------------------------------------------

  for r = 1:nstart
      fprintf('-----------------------------------------------------------------------------------------\n');
      fprintf('         (b.%d)   NPL Algorithm: %s\n', r, start_name{r});
      [ theta_npl{r}, varb, pnpl, start_v1npl{r}, vcnpl, start_ll(r), start_iter(r), start_err(r), start_times{r} ] = ...
          npldygam1(aobs, zexo_obs, aobs_1, z_exo_val, ptranz_exo_est, start_p{r}, disfact_est, miniter, maxiter_npl);

       start_theta1(:,r)=reshape(theta_npl{r}(1,:),[],1);%%%%%
        %%start_theta1(:,r)=best_logit;
  end

  max_llike = max(start_ll);
  best_r = find(start_ll == max_llike, 1, 'first');
  if best_r && sum(start_err) < nstart
      bmat_1npl(draw,:) = theta_npl{best_r}(1,:);
      bmat_2npl(draw,:) = theta_npl{best_r}(2,:);
      bmat_3npl(draw,:) = theta_npl{best_r}(3,:);
      bmat_cnpl(draw,:) = theta_npl{best_r}(start_iter(best_r),:);
      time_1npl(draw,1) = start_times{best_r}(1);
      time_2npl(draw,1) = sum(start_times{best_r}(1:2));
      time_3npl(draw,1) = sum(start_times{best_r}(1:3));
      time_cnpl(draw,1) = sum(start_times{best_r});
      iter_cnpl(draw,1) = start_iter(best_r);
      time_2npl_iter(draw,1) = time_2npl(draw,1) / 2;
      time_3npl_iter(draw,1) = time_3npl(draw,1) / 3;
      time_cnpl_iter(draw,1) = time_cnpl(draw,1) / iter_cnpl(draw,1);
      fail_cnpl(draw,1) = 0;
  else
      fail_cnpl(draw,1) = 1;
  end


  %--------------------------------------------------
  % EPL
  %--------------------------------------------------

    for r = 1:nstart
      fprintf('-----------------------------------------------------------------------------------------\n');
      fprintf('         (b.%d)   EPL Algorithm: %s\n', r, start_name{r});
      if (start_1npl(r) > 0)
          v0 = start_v1npl{r};
          theta0 = start_theta1(:,r);
      else
          v0 = start_v{r};
          theta0 = start_theta(:,r);
      end
      

      %% EPL Using analytical Jacobian

      ww_array=[];

      %compute_jacobian_spec=1;%%%%%
      %krylov_spec=0;%%%%%

      tStart = tic; 
      [ theta_epl{r}, varb, start_ll(r), start_iter(r), start_err(r), start_times{r},ww_array{r},vsol_epl ] =...
          epldygam(aobs, zexo_obs, aobs_1, z_exo_val, ptranz_exo_est, theta0, v0, disfact_est, miniter, maxiter);
      toc(tStart)

      if isempty(ww_array)==0
          ww_array_stack(:,:,:,draw)=ww_array{1};
      end
      param_epl=reshape(theta_epl{r}(end,:),[],1);
    end

    max_llike = max(start_ll);
    best_r = find(start_ll == max_llike, 1, 'first');
    if best_r && sum(start_err) < nstart
        bmat_1epl(draw,:) = theta_epl{best_r}(1,:);
        bmat_2epl(draw,:) = theta_epl{best_r}(2,:);
        bmat_3epl(draw,:) = theta_epl{best_r}(3,:);
        bmat_cepl(draw,:) = theta_epl{best_r}(start_iter(best_r),:);
        time_1epl(draw,1) = start_times{best_r}(1);
        time_2epl(draw,1) = sum(start_times{best_r}(1:2));
        time_3epl(draw,1) = sum(start_times{best_r}(1:3));
        time_cepl(draw,1) = sum(start_times{best_r});
        iter_cepl(draw,1) = start_iter(best_r);
        time_2epl_iter(draw,1) = time_2epl(draw,1) / 2;
        time_3epl_iter(draw,1) = time_3epl(draw,1) / 3;
        time_cepl_iter(draw,1) = time_cepl(draw,1) / iter_cepl(draw,1);
        fail_cepl(draw,1) = 0;
    else
        fail_cepl(draw,1) = 1;
    end
  
    if 1==1 % sequential
        run run_sequential_dynamic_game.m
    end

end% draw
    


% True Parameters
theta_true = trueparam(1:kparam,1)';

% Replication selection
keep_cnpl = fail_cnpl == 0;
keep_cepl = fail_cepl == 0;

% Empirical Means
mean_1npl = mean(bmat_1npl);
mean_2npl = mean(bmat_2npl);
mean_3npl = mean(bmat_3npl);
mean_cnpl = mean(bmat_cnpl(keep_cnpl,:));
mean_1epl = mean(bmat_1epl);
mean_2epl = mean(bmat_2epl);
mean_3epl = mean(bmat_3epl);
mean_cepl = mean(bmat_cepl(keep_cepl,:));
mean_1npl_true = mean(btrue_npl);
mean_1epl_true = mean(btrue_epl);

% Mean Bias
meanbias_1npl = mean_1npl - theta_true;
meanbias_2npl = mean_2npl - theta_true;
meanbias_3npl = mean_3npl - theta_true;
meanbias_cnpl = mean_cnpl - theta_true;
meanbias_1epl = mean_1epl - theta_true;
meanbias_2epl = mean_2epl - theta_true;
meanbias_3epl = mean_3epl - theta_true;
meanbias_cepl = mean_cepl - theta_true;
meanbias_1npl_true = mean_1npl_true - theta_true;
meanbias_1epl_true = mean_1epl_true - theta_true;

% Empirical Variances
var_1npl = var(bmat_1npl);
var_2npl = var(bmat_2npl);
var_3npl = var(bmat_3npl);
var_cnpl = var(bmat_cnpl(keep_cnpl,:));
var_1epl = var(bmat_1epl);
var_2epl = var(bmat_2epl);
var_3epl = var(bmat_3epl);
var_cepl = var(bmat_cepl(keep_cepl,:));
var_1npl_true = var(btrue_npl);
var_1epl_true = var(btrue_epl);

% MSE
mse_1npl = meanbias_1npl.^2 + var_1npl;
mse_2npl = meanbias_2npl.^2 + var_2npl;
mse_3npl = meanbias_3npl.^2 + var_3npl;
mse_cnpl = meanbias_cnpl.^2 + var_cnpl;
mse_1epl = meanbias_1epl.^2 + var_1epl;
mse_2epl = meanbias_2epl.^2 + var_2epl;
mse_3epl = meanbias_3epl.^2 + var_3epl;
mse_cepl = meanbias_cepl.^2 + var_cepl;
mse_1npl_true = meanbias_1npl_true.^2 + var_1npl_true;
mse_1epl_true = meanbias_1epl_true.^2 + var_1epl_true;

% Iterations
median_iter_cnpl = median(iter_cnpl(keep_cnpl));
median_iter_cepl = median(iter_cepl(keep_cepl));
max_iter_cnpl = max(iter_cnpl(keep_cnpl));
max_iter_cepl = max(iter_cepl(keep_cepl));
iqr_iter_cnpl = iqr(iter_cnpl(keep_cnpl));
iqr_iter_cepl = iqr(iter_cepl(keep_cepl));
iter_fail_cnpl = sum(iter_cnpl(keep_cnpl) == maxiter);
iter_fail_cepl = sum(iter_cepl(keep_cepl) == maxiter);

% Total time and time per iteration
median_time_1npl = median(time_1npl(keep_cnpl));
median_time_1epl = median(time_1epl(keep_cepl));
median_time_2npl = median(time_2npl(keep_cnpl));
median_time_2epl = median(time_2epl(keep_cepl));
median_time_3npl = median(time_3npl(keep_cnpl));
median_time_3epl = median(time_3epl(keep_cepl));
median_time_cnpl = median(time_cnpl(keep_cnpl));
median_time_cepl = median(time_cepl(keep_cepl));

median_time_1npl_iter = median_time_1npl;
median_time_1epl_iter = median_time_1epl;
median_time_2npl_iter = median(time_2npl_iter(keep_cnpl));
median_time_2epl_iter = median(time_2epl_iter(keep_cepl));
median_time_3npl_iter = median(time_3npl_iter(keep_cnpl));
median_time_3epl_iter = median(time_3epl_iter(keep_cepl));
median_time_cnpl_iter = median(time_cnpl_iter(keep_cnpl));
median_time_cepl_iter = median(time_cepl_iter(keep_cepl));

mean_time_1npl = mean(time_1npl(keep_cnpl));
mean_time_1epl = mean(time_1epl(keep_cepl));
mean_time_2npl = mean(time_2npl(keep_cnpl));
mean_time_2epl = mean(time_2epl(keep_cepl));
mean_time_3npl = mean(time_3npl(keep_cnpl));
mean_time_3epl = mean(time_3epl(keep_cepl));
mean_time_cnpl = mean(time_cnpl(keep_cnpl));
mean_time_cepl = mean(time_cepl(keep_cepl));

mean_time_1npl_iter = mean_time_1npl;
mean_time_1epl_iter = mean_time_1epl;
mean_time_2npl_iter = mean(time_2npl_iter(keep_cnpl));
mean_time_2epl_iter = mean(time_2epl_iter(keep_cepl));
mean_time_3npl_iter = mean(time_3npl_iter(keep_cnpl));
mean_time_3epl_iter = mean(time_3epl_iter(keep_cepl));
mean_time_cnpl_iter = mean(time_cnpl_iter(keep_cnpl));
mean_time_cepl_iter = mean(time_cepl_iter(keep_cepl));

total_time_1npl = sum(time_1npl(keep_cnpl));
total_time_1epl = sum(time_1epl(keep_cepl));
total_time_2npl = sum(time_2npl(keep_cnpl));
total_time_2epl = sum(time_2epl(keep_cepl));
total_time_3npl = sum(time_3npl(keep_cnpl));
total_time_3epl = sum(time_3epl(keep_cepl));
total_time_cnpl = sum(time_cnpl(keep_cnpl));
total_time_cepl = sum(time_cepl(keep_cepl));

table1=[theta_true', ...
            meanbias_1npl', meanbias_1epl', meanbias_2npl', meanbias_2epl', ...
            meanbias_3npl', meanbias_3epl', meanbias_cnpl', meanbias_cepl'];
table2=[theta_true', ...
            mse_1npl', mse_1epl', mse_2npl', mse_2epl', ...
            mse_3npl', mse_3epl', mse_cnpl', mse_cepl'];
table3=[[NaN;NaN;NaN;NaN], ...
            [1;1;0;0], [1;1;0;0], [2;2;0;0], [2;2;0;0], ...
            [3;3;0;0], [3;3;0;0], ...
            [median_iter_cnpl;max_iter_cnpl;iqr_iter_cnpl;iter_fail_cnpl],...
            [median_iter_cepl;max_iter_cepl;iqr_iter_cepl;iter_fail_cepl]];
table4=[[NaN;NaN;NaN;NaN],...
    [total_time_1npl;mean_time_1npl;median_time_1npl;mean_time_1npl_iter],...
    [total_time_1epl;mean_time_1epl;median_time_1epl;mean_time_1epl_iter],...
    [total_time_2npl;mean_time_2npl;median_time_2npl;mean_time_2npl_iter],...
    [total_time_2epl;mean_time_2epl;median_time_2epl;mean_time_2epl_iter],...
    [total_time_3npl;mean_time_3npl;median_time_3npl;mean_time_3npl_iter],...
    [total_time_3epl;mean_time_3epl;median_time_3epl;mean_time_3epl_iter],...
    [total_time_cnpl;mean_time_cnpl;median_time_cnpl;mean_time_cnpl_iter],...
    [total_time_cepl;mean_time_cepl;median_time_cepl;mean_time_cepl_iter]];
    
table=[table1;table2;table3;table4];



