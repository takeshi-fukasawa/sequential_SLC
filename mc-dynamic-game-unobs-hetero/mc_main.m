% mc_main.m -- Control program for running Monte Carlo experiments.


clear all
global selexper
global nobs
global param_seq iter_info_sequential
global compute_jacobian_spec krylov_spec
global flag_vec relres_cell iter_gmres_cell resvec_cell
global output_opt
global rho_val_global 
global VFI_spec

warning('off', 'all')

sequential_code_path='C:/Users/fukas/Dropbox/git/sequential/';
run(append(sequential_code_path,'addpath_preliminary.m'))

addpath('./results')
addpath('./unobs_hetero')

save_spec=1;

selexper=3;
nobs=6400;
%nobs=1600;
unobs_hetero_spec=1;
analytical_spec=1;% Use analytical derivatives
nrepli = 20;%20;%1000  % Number of Monte carlo replications
nplayer = 3;    %  Number of players
n_sim=nrepli;
nstart_seq=5;
VFI_spec=0;
spectral_spec=1;
run_nfxp_spec=0;

if unobs_hetero_spec==0
    multiple_zexo_spec=0;
else
    multiple_zexo_spec=1;
end


compute_jacobian_spec=1;
krylov_spec=0;

%%%%%%%%%%%
obj_val_SLC_all=NaN(n_sim,nstart_seq);
times_SLC_all=NaN(n_sim,nstart_seq);
n_iter_SLC_all=NaN(n_sim,nstart_seq);
count_G_eval_SLC_all=NaN(n_sim,nstart_seq);
count_obj_eval_SLC_all=NaN(n_sim,nstart_seq);

obj_val_nfxp_all=NaN(n_sim,nstart_seq);
times_nfxp_all=NaN(n_sim,nstart_seq);
n_iter_nfxp_all=NaN(n_sim,nstart_seq);
count_G_eval_nfxp_all=NaN(n_sim,nstart_seq);
count_obj_eval_nfxp_all=NaN(n_sim,nstart_seq);
%%%%%%%%%%%%%

mc_sequential
