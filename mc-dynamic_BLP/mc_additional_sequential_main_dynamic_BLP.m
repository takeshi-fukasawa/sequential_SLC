clear all
code_path='C:/Users/fukas/Dropbox/git/';
run(append(code_path,'sequential_SLC_temp/addpath_preliminary.m'))

run preliminary_set_params_dynamic_BLP.m

check_conv_to_same_param_spec=0;
run_nfxp_spec=0;
run_EPL_spec=0;
spectral_spec=0;
SLC_tol=1e-6;


save_spec=1;
n_market=20;
nstart=5;

%% SLC with bad initial values
n_market=1;
nstart=1;
bad_init_val_spec=0;
spectral_spec=0;
SLC_line_search_spec=1;
%run run_sequential_dynamic_BLP.m

%% SLC with line search
check_conv_to_same_param_spec=0;
bad_init_val_spec=0;
SLC_line_search_spec=1;
spectral_spec=0;
n_market=20;
nstart=5;

%run run_sequential_dynamic_BLP.m

%% SLC-spectral with line search
check_conv_to_same_param_spec=0;
bad_init_val_spec=0;
SLC_line_search_spec=1;
spectral_spec=1;
n_market=20;
nstart=5;

%run run_sequential_dynamic_BLP.m

%% SLC-spectral with more initial values
check_conv_to_same_param_spec=0;
bad_init_val_spec=0;
SLC_line_search_spec=0;
spectral_spec=1;
n_market=20;
nstart=10;

run run_sequential_dynamic_BLP.m
