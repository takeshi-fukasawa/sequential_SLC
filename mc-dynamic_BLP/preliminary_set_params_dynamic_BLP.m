BLP_algorithm_code_path=append(code_path,'BLP_algorithm/');

addpath(append(BLP_algorithm_code_path,'static_dynamic_BLP_algorithm'))


%clear

global Pr0_spec inv_multiply_Chebyshev_basis Chebyshev_extrema
global ratio_s_i0t_ccp_t
global V_initial_global
global positive_param_nfxp_spec

rng('default')


addpath(append(BLP_algorithm_code_path,'static_dynamic_BLP_algorithm/functions'))
addpath(append(BLP_algorithm_code_path,'static_dynamic_BLP_algorithm/run_dynamic'))
addpath(append(BLP_algorithm_code_path,'static_dynamic_BLP_algorithm/run_static_dynamic'))


spec_default=[];
spec_default.spectral_coef_spec=3;
spec_default.TOL=1e-12;
spec_default.common_spectral_coef_spec=0;
spec_default.ITER_MAX=3000;
%spec_default.ITER_MAX=3;
%%spec_default.stabilize_spec=1;
spec_default.DEBUG=1;

skip_contraction_spec=0;
large_hetero_spec=0;
mistake_spec=1;%%%%

rho_true=0;
rho_est=rho_true;
L=100000;

G=1;

n_sim=1;
I=50;
ns=I;

n_draw=1;
initial_scale=2;

IVS_spec=1; % If IVS_spec==1, use inclusive value sufficiency. If IVS_spec==0, assume perfect foresight on the future state transitions of market-level variables. For details, see the paper.

%%%%%%%%%%%%%%%%
%% Simulation 1
durable_spec=1; % Use the model in Sun and Ishihara (2019)
%%durable_spec=0;

Pr0_spec=1; %%%% Pr0_spec==1=> Introduce endogenous Pr0

t_dependent_spectral_coef_spec=0; % Introduce time-dependent step sizes alpha in the spectral / SQUAREM algorithms.

n_dim_V=1;

J=25;% Number of products per nest

beta_C=0.99;
%beta_C=0.01;

%%% T==1=> stationary;

if IVS_spec==1
   T=25;
   n_grid_IV=50;
   %n_grid_IV=20;
   beta_0=6;
   Pr0_spec=1;
   t_dependent_spectral_coef_spec=0;%%%
   %spec_default.spectral_coef_max=5;%%%%
   spec_default.ITER_MAX=3000;
   n_draw=5;
else
    T=50;
    n_grid_IV=0;
    beta_0=6;
    Pr0_spec=1;
    t_dependent_spectral_coef_spec=1;
    spec_default.spectral_coef_max=10;
    spec_default.dim_hetero_alpha=4;
end

%%T=1;


n_col_results=5;

TOL_DIST_s_jt=1e-12;
t_dim_id=4;

n_dim_V=T+n_grid_IV;
tune_param=1;Newton_spec=0;

%%%%%%%%%%%%%%%%%%
%% Parameter settings for sequential estimation testing
n_market=1;
n_sim=n_market;
market_id=1;

nstart=1;

file_name_preliminary_set_params='preliminary_set_params_dynamic_BLP.m';

positive_param_nfxp_spec=0;
