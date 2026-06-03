%% Run NFXP, SLC
check_conv_to_same_param_spec=0;
run_nfxp_spec=1;
spectral_spec=0;
SLC_tol=1e-6;
run run_sequential_dynamic_BLP.m

%% Run SLC(loose tolerance)
check_conv_to_same_param_spec=0;
run_nfxp_spec=0;
spectral_spec=0;
SLC_tol=1e-4;
run run_sequential_dynamic_BLP.m

%% Run SLC-spectral
check_conv_to_same_param_spec=0;
run_nfxp_spec=0;
spectral_spec=1;
SLC_tol=1e-6;
run run_sequential_dynamic_BLP.m

if 1==0
    %% Check convergence (SLC; from perturbed NFXP estimate)
    check_conv_to_same_param_spec=1;
    run_nfxp_spec=0;
    spectral_spec=0;
    SLC_tol=1e-6;
    run run_sequential_dynamic_BLP.m

    ratio_SLC_NFXP=reshape(param_SLC_all,nrepli,n_param)./param_nfxp_est
    
end
