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

code_path='C:/Users/fukas/Dropbox/git/';
run(append(code_path,'sequential_SLC_temp/addpath_preliminary.m'))



addpath(append(sequential_code_path,'mc-dynamic-game-unobs-hetero/results'))
addpath(append(sequential_code_path,'mc-dynamic-game-unobs-hetero/unobs_hetero'))


save_spec=0;

selexper=3;
unobs_hetero_spec=1;
nrepli = 20;%20;%1000  % Number of Monte carlo replications
nplayer = 3;    %  Number of players
n_sim=nrepli;
nstart_seq=5;
spectral_spec=0;
run_nfxp_spec=0;

if unobs_hetero_spec==0
    multiple_zexo_spec=0;
else
    multiple_zexo_spec=1;
end

compute_jacobian_spec=1;
krylov_spec=0;

results_table_SLC_EPL_comparison=[];

for mc_additional_type_id=3:3
    if mc_additional_type_id==1
        % SLC vs EPL (Analytical derivativesa)
        nobs=6400;
        VFI_spec=0;
        nrepli=1;nstart_seq=1;
        run_NFXP_spec=0;
        run_EPL_spec=1;
        analytical_spec=1;
        SLC_line_search_spec=0;
        bad_init_val_spec=0;
    elseif mc_additional_type_id==2
       % SLC vs EPL (Numerical derivatives)
        nobs=6400;
        VFI_spec=0;
        nrepli=1;nstart_seq=1;
        run_NFXP_spec=0;
        run_EPL_spec=1;
        analytical_spec=0;
        SLC_line_search_spec=0;
        bad_init_val_spec=0;
    elseif mc_additional_type_id==3
     % Bad initial values (Analytical)
        nobs=400;
        VFI_spec=0;
        nrepli=1;nstart_seq=1;
        run_NFXP_spec=0;
        run_EPL_spec=0;
        analytical_spec=1;
        SLC_line_search_spec=0;
        bad_init_val_spec=0;
    end

     n_sim=nrepli;

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

    run mc_sequential.m

    if mc_additional_type_id<=2
        iter_info_EPL=results_EPL.iter_info_EPL;
    
        temp_SLC=[iter_info_SLC.t_cpu,iter_info_SLC.count_G_eval,...
            iter_info_SLC.count_obj_eval,iter_info_SLC.n_iter];
        temp_EPL=[iter_info_EPL.t_cpu,iter_info_EPL.count_G_eval,...
            iter_info_EPL.count_obj_eval,iter_info_EPL.n_iter];
        
        results_table_SLC_EPL_comparison=[results_table_SLC_EPL_comparison;...
            temp_SLC;temp_EPL];
    end
    
    if save_spec==1 & mc_additional_type_id==2
        file_name=append(save_path,"dynamic_game_unobs_hetero_SLC_EPL_comparison.csv");
        writematrix(round(results_table_SLC_EPL_comparison,1),file_name);   
    end

end
