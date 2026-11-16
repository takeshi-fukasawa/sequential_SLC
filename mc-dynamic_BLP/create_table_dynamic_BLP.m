code_path='C:/Users/fukas/Dropbox/git/';
run(append(code_path,'sequential_SLC_temp/addpath_preliminary.m'))

run preliminary_set_params_dynamic_BLP.m

nrepli=20;
sim_id_max=1;

n_param=3;


load(append(save_path,"dynamic_BLP_seq_results.mat"));%,...
%%% Including NFXP results 

[obj_val_SLC_all,param_SLC_all,iter_info_SLC_all,...
    times_SLC_all,...
    count_G_eval_SLC_all,count_obj_eval_SLC_all,iter_SLC_all,...
    feval_SLC_all,div_SLC_all]=...
    read_saved_data_func("dynamic_BLP_seq_results.mat",save_path);

[obj_val_SLC_spectral_all,param_SLC_spectral_all,iter_info_SLC_spectral_all,...
    times_SLC_spectral_all,...
    count_G_eval_SLC_spectral_all,count_obj_eval_SLC_spectral_all,iter_SLC_spectral_all,...
    feval_SLC_spectral_all,div_SLC_spectral_all]=...
    read_saved_data_func("dynamic_BLP_seq_results_SLC_only_spectral.mat",save_path);

[obj_val_SLC_spectral_bad_init_all,param_SLC_spectral_bad_init_all,iter_info_SLC_spectral_bad_init_all,...
    times_SLC_spectral_bad_init_all,...
    count_G_eval_SLC_spectral_bad_init_all,count_obj_eval_SLC_spectral_bad_init_all,iter_SLC_spectral_bad_init_all,...
    feval_SLC_spectral_bad_init_all,div_SLC_spectral_bad_init_all]=...
    read_saved_data_func("dynamic_BLP_seq_results_SLC_only_spectral_bad_init.mat",save_path);

[obj_val_SLC_spectral_many_starting_values_all,param_SLC_spectral_many_starting_values_all,iter_info_SLC_spectral_many_starting_values_all,...
    times_SLC_spectral_many_starting_values_all,...
    count_G_eval_SLC_spectral_many_starting_values_all,count_obj_eval_SLC_spectral_many_starting_values_all,iter_SLC_spectral_many_starting_values_all,...
    feval_SLC_spectral_many_starting_values_all,div_SLC_spectral_many_starting_values_all]=...
    read_saved_data_func("dynamic_BLP_seq_results_SLC_only_spectral_many_starting_values.mat",save_path);


[obj_val_SLC_loose_tol_all,param_SLC_loose_tol_all,iter_info_SLC_loose_tol_all,...
    times_SLC_loose_tol_all,...
    count_G_eval_SLC_loose_tol_all,count_obj_eval_SLC_loose_tol_all,iter_SLC_loose_tol_all,...
    feval_SLC_loose_tol_all,div_SLC_loose_tol_all]=...
    read_saved_data_func("dynamic_BLP_seq_results_SLC_only_SLC_loose_tol.mat",save_path);


[obj_val_SLC_line_search_all,param_SLC_line_search_all,iter_info_SLC_line_search_all,...
    times_SLC_line_search_all,...
    count_G_eval_SLC_line_search_all,count_obj_eval_SLC_line_search_all,iter_SLC_line_search_all,...
    feval_SLC_line_search_all,div_SLC_line_search_all]=...
    read_saved_data_func("dynamic_BLP_seq_results_SLC_only_line_search.mat",save_path);

[obj_val_SLC_spectral_line_search_all,param_SLC_spectral_line_search_all,iter_info_SLC_spectral_line_search_all,...
    times_SLC_spectral_line_search_all,...
    count_G_eval_SLC_spectral_line_search_all,count_obj_eval_SLC_spectral_line_search_all,iter_SLC_spectral_line_search_all,...
    feval_SLC_spectral_line_search_all,div_SLC_spectral_line_search_all]=...
    read_saved_data_func("dynamic_BLP_seq_results_SLC_only_spectral_line_search.mat",save_path);

for sim_id=1:sim_id_max


    if sim_id==1
        results_param_est=cell(sim_id_max,1);
        results_comp=cell(sim_id_max,1);

        param_nfxp_est_array=NaN(nrepli,n_param,3);
        param_SLC_est_array=NaN(nrepli,n_param,3);
        param_SLC_spectral_est_array=NaN(nrepli,n_param,3);
        
    end

    MAX_ITER_SLC=50;
    non_conv_SLC=(iter_SLC_all+1>=MAX_ITER_SLC|iter_SLC_all==1|div_SLC_all==1);
    non_conv_SLC_spectral=(iter_SLC_spectral_all+1>=MAX_ITER_SLC|iter_SLC_spectral_all==1|div_SLC_spectral_all==1);
    non_conv_SLC_loose_tol=(iter_SLC_loose_tol_all+1>=MAX_ITER_SLC|iter_SLC_loose_tol_all==1|div_SLC_loose_tol_all==1);
    non_conv_SLC_line_search=(iter_SLC_line_search_all+1>=MAX_ITER_SLC|iter_SLC_line_search_all==1|div_SLC_line_search_all==1);
    non_conv_SLC_spectral_line_search=(iter_SLC_spectral_line_search_all+1>=MAX_ITER_SLC|iter_SLC_spectral_line_search_all==1|div_SLC_spectral_line_search_all==1);
    non_conv_SLC_spectral_bad_init=(iter_SLC_spectral_bad_init_all+1>=MAX_ITER_SLC|iter_SLC_spectral_bad_init_all==1|div_SLC_spectral_bad_init_all==1);
    non_conv_SLC_spectral_many_starting_values=(iter_SLC_spectral_many_starting_values_all+1>=MAX_ITER_SLC|iter_SLC_spectral_many_starting_values_all==1|div_SLC_spectral_many_starting_values_all==1);
    MAX_ITER_nfxp=100;%%%
    non_conv_nfxp=(iter_nfxp_all>=MAX_ITER_nfxp|iter_nfxp_all==1);
    obj_val_SLC_all(non_conv_SLC)=1e10;
    obj_val_SLC_spectral_all(non_conv_SLC_spectral)=1e10;
    obj_val_SLC_loose_tol_all(non_conv_SLC_loose)=1e10;
    obj_val_SLC_line_search_all(non_conv_SLC_line_search)=1e10;
    obj_val_SLC_spectral_many_starting_values_all(non_conv_SLC_spectral_many_starting_values)=1e10;
    obj_val_SLC_spectral_line_search_all(non_conv_SLC_spectral_line_search)=1e10;
    obj_val_nfxp_all(non_conv_nfxp)=1e10;


    [param_SLC_est,min_val_SLC]=obtain_param_SLC_func(param_SLC_all,obj_val_SLC_all);

    [param_SLC_spectral_est,min_val_SLC_spectral]=obtain_param_SLC_func(param_SLC_spectral_all,obj_val_SLC_spectral_all);
    [param_SLC_spectral_many_starting_values_est,min_val_SLC_spectral_many_starting_values]=obtain_param_SLC_func(param_SLC_spectral_many_starting_values_all,obj_val_SLC_spectral_many_starting_values_all);

    [param_SLC_loose_tol_est,min_val_SLC_loose_tol]=obtain_param_SLC_func(param_SLC_loose_tol_all,obj_val_SLC_loose_tol_all);
    [param_SLC_line_search_est,min_val_SLC_line_search]=obtain_param_SLC_func(param_SLC_line_search_all,obj_val_SLC_line_search_all);

    [param_SLC_spectral_line_search_est,min_val_SLC_spectral_line_search]=obtain_param_SLC_func(param_SLC_spectral_line_search_all,obj_val_SLC_spectral_line_search_all);

    [param_SLC_spectral_bad_init_est,min_val_SLC_spectral_bad_init]=obtain_param_SLC_func(param_SLC_spectral_bad_init_all,obj_val_SLC_spectral_bad_init_all);

    [param_nfxp_est,min_val_nfxp]=obtain_param_SLC_func(param_nfxp_all,obj_val_nfxp_all);

    
    times_nfxp_sum=sum(times_nfxp_all,2);
    times_SLC_sum=sum(times_SLC_all,2);
    times_SLC_spectral_sum=sum(times_SLC_spectral_all,2);
    times_SLC_loose_tol_sum=sum(times_SLC_loose_tol_all,2);
    times_SLC_line_search_sum=sum(times_SLC_line_search_all,2);
    times_SLC_spectral_line_search_sum=sum(times_SLC_spectral_line_search_all,2);
    times_SLC_spectral_bad_init_sum=sum(times_SLC_spectral_bad_init_all,2);
    %%%%
    

    col_nfxp=construct_table_col_func(times_nfxp_all,count_obj_eval_nfxp_all,count_G_eval_nfxp_all,iter_nfxp_all,non_conv_nfxp);

    col_SLC=construct_table_col_func(times_SLC_all,count_obj_eval_SLC_all,count_G_eval_SLC_all,iter_SLC_all,non_conv_SLC);
    col_SLC_spectral=construct_table_col_func(times_SLC_spectral_all,count_obj_eval_SLC_spectral_all,count_G_eval_SLC_spectral_all,iter_SLC_spectral_all,non_conv_SLC_spectral);
    col_SLC_loose=construct_table_col_func(times_SLC_loose_tol_all,count_obj_eval_SLC_loose_tol_all,count_G_eval_SLC_loose_tol_all,iter_SLC_loose_tol_all,non_conv_SLC_loose_tol);
    col_SLC_line_search=construct_table_col_func(times_SLC_line_search_all,count_obj_eval_SLC_line_search_all,count_G_eval_SLC_line_search_all,iter_SLC_line_search_all,non_conv_SLC_line_search);
    col_SLC_spectral_line_search=construct_table_col_func(times_SLC_spectral_line_search_all,count_obj_eval_SLC_spectral_line_search_all,count_G_eval_SLC_spectral_line_search_all,iter_SLC_spectral_line_search_all,non_conv_SLC_spectral_line_search);
    col_SLC_spectral_bad_init=construct_table_col_func(times_SLC_spectral_bad_init_all,count_obj_eval_SLC_spectral_bad_init_all,count_G_eval_SLC_spectral_bad_init_all,iter_SLC_spectral_bad_init_all,non_conv_SLC_spectral_bad_init);
    col_SLC_spectral_many_starting_values=construct_table_col_func(times_SLC_spectral_many_starting_values_all,count_obj_eval_SLC_spectral_many_starting_values_all,count_G_eval_SLC_spectral_many_starting_values_all,iter_SLC_spectral_many_starting_values_all,non_conv_SLC_spectral_many_starting_values);


    results_comp_table=[col_SLC;col_SLC_loose;col_SLC_spectral;col_nfxp];
    results_param_est=[mean(abs(param_SLC_est));std(abs(param_SLC_est));...
            mean(abs(param_SLC_spectral_est));std(abs(param_SLC_spectral_est));...
        mean(abs(param_nfxp_est));std(abs(param_nfxp_est))];

    param_nfxp_est_array(:,:,sim_id)=param_nfxp_est;
    param_SLC_est_array(:,:,sim_id)=param_SLC_est;
    param_SLC_loose_tol_est_array(:,:,sim_id)=param_SLC_loose_tol_est;
    param_SLC_spectral_est_array(:,:,sim_id)=param_SLC_spectral_est;
    param_SLC_line_search_est_array(:,:,sim_id)=param_SLC_line_search_est;
    param_SLC_spectral_line_search_est_array(:,:,sim_id)=param_SLC_spectral_line_search_est;
    param_SLC_spectral_bad_init_est_array(:,:,sim_id)=param_SLC_spectral_bad_init_est;
    param_SLC_spectral_many_starting_values_est_array(:,:,sim_id)=param_SLC_spectral_many_starting_values_est;
end



writematrix(round(results_comp_table,1),append(save_path,"dynamic_BLP_table_comp.csv"))
writematrix(round(results_param_est,3),append(save_path,"dynamic_BLP_table_param_est.csv"))

min_val_difference_SLC=min(min_val_SLC,[],2)-min(min_val_nfxp,[],2)
min_val_difference_SLC_spectral=min(min_val_SLC_spectral,[],2)-min(min_val_nfxp,[],2)
min_val_difference_SLC_spectral_many_starting_values=min(min_val_SLC_spectral_many_starting_values,[],2)-min(min_val_nfxp,[],2)

diff_est_NFXP_SLC_loose=abs(param_SLC_loose_tol_est)-abs(param_nfxp_est);
diff_est_NFXP_SLC_spectral=abs(param_SLC_spectral_est)-abs(param_nfxp_est);
diff_est_NFXP_SLC_spectral_many_starting_values=abs(param_SLC_spectral_many_starting_values_est)-abs(param_nfxp_est);

mean_feval_SLC_line_search=mean(feval_SLC_line_search_all(:)+1);
mean_feval_SLC_spectral_line_search=mean(feval_SLC_spectral_line_search_all(:)+1);

results_comp_table_line_search=[col_SLC;col_SLC_spectral;...
    col_SLC_line_search;col_SLC_spectral_line_search];
writematrix(round(results_comp_table_line_search,1),...
    append(save_path,"dynamic_BLP_table_comp_line_search.csv"))

mean_feval_line_search_vec=[mean_feval_SLC_line_search;mean_feval_SLC_spectral_line_search]

results_comp_table_bad_init=[col_SLC_spectral;col_SLC_spectral_bad_init];
results_comp_table_bad_init(1,1:2)=results_comp_table_bad_init(1,1:2)./nstart;% Show average time per starting value for the first row (SLC_spectral) to be consistent with the second row (SLC_spectral_bad_init) in terms of number of starting values.
writematrix(round(results_comp_table_bad_init,1),...
    append(save_path,"dynamic_BLP_table_comp_bad_init.csv"))
