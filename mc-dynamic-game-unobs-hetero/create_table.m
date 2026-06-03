run("C:/Users/fukas/Dropbox/git/sequential/addpath_preliminary.m")
save_path="C:/Users/fukas/Dropbox/sequential/output/";


sim_id_max=5;
for sim_id=1:sim_id_max
    if sim_id==1
        file_name=append(save_path,"dynamic_game_unobs_hetero_seq_results_6400_analytical_diff.mat");
    elseif sim_id==2
        file_name=append(save_path,"dynamic_game_unobs_hetero_seq_results_numerical_diff.mat");
    elseif sim_id==3
        file_name=append(save_path,"dynamic_game_unobs_hetero_seq_results_1600_analytical_diff.mat");
    elseif sim_id==4
        file_name=append(save_path,"dynamic_game_unobs_hetero_seq_results_analytical_diff_EVFI.mat");
    elseif sim_id==5
        file_name=append(save_path,"dynamic_game_unobs_hetero_seq_results_numerical_diff_EVFI.mat");
    end

    load(file_name)

    if sim_id==1
        results_param_est=cell(sim_id_max,1);
        results_comp=cell(sim_id_max,1);

        param_nfxp_est_array=NaN(nrepli,n_param,3);
        param_SLC_est_array=NaN(nrepli,n_param,3);
    end

    MAX_ITER_SLC=50;
    non_conv_SLC=(iter_SLC_all+1>=MAX_ITER_SLC);
    MAX_ITER_nfxp=100;%%%
    non_conv_nfxp=(iter_nfxp_all>=MAX_ITER_nfxp);
    obj_val_SLC_all(non_conv_SLC)=1e10;
    obj_val_nfxp_all(non_conv_nfxp)=1e10;

    [min_val, temp] = min(obj_val_SLC_all, [], 2);
    param_SLC_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_SLC_all_i = param_SLC_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_SLC_all_i), rows, temp);
        param_SLC_est(:,i) = param_SLC_all_i(idx);
    end
    
    [min_val, temp] = min(obj_val_nfxp_all, [], 2);
    param_nfxp_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_nfxp_all_i = param_nfxp_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_nfxp_all_i), rows, temp);
        param_nfxp_est(:,i) = param_nfxp_all_i(idx);
    end
    
    times_nfxp_sum=sum(times_nfxp_all,2);
    times_SLC_sum=sum(times_SLC_all,2);
    
    %%%%
    speed_test_info
    
    col_nfxp=[results_time_dist_func(sum(times_nfxp_all,2)),...
        mean(count_obj_eval_nfxp_all(:),'omitnan'),...
        mean(count_G_eval_nfxp_all(:),'omitnan'),...
        mean(iter_nfxp_all(:)),...
        mean(1-min(non_conv_nfxp,[],2))*100,mean(1-non_conv_nfxp(:))*100];
    col_SLC=[results_time_dist_func(sum(times_SLC_all,2)),...
        mean(count_obj_eval_SLC_all(:),'omitnan'),...
        mean(count_G_eval_SLC_all(:),'omitnan'),...
        mean(iter_SLC_all(:)+1),...
        mean(1-min(non_conv_SLC,[],2))*100,mean(1-non_conv_SLC(:))*100];
    results_comp{sim_id}=[col_SLC;col_nfxp];
    results_param_est{sim_id}=[mean(param_SLC_est);std(param_SLC_est);...
        mean(param_nfxp_est);std(param_nfxp_est)];

    param_nfxp_est_array(:,:,sim_id)=param_nfxp_est;
    param_SLC_est_array(:,:,sim_id)=param_SLC_est;
    
end

param_diff_6400_nfxp_analy=max(abs(param_nfxp_est_array(:,:,1)-param_SLC_est_array(:,:,1)),[],'all');
param_diff_6400_SLC_numer=max(abs(param_SLC_est_array(:,:,2)-param_SLC_est_array(:,:,1)),[],'all');
param_diff_6400_nfxp_numer=max(abs(param_nfxp_est_array(:,:,2)-param_SLC_est_array(:,:,1)),[],'all');
param_diff_6400_SLC_numer_EVFI=max(abs(param_SLC_est_array(:,:,5)-param_SLC_est_array(:,:,1)),[],'all');
param_diff_6400_SLC_analy_EVFI=max(abs(param_SLC_est_array(:,:,4)-param_SLC_est_array(:,:,1)),[],'all');
param_diff_6400_nfxp_analy_EVFI=max(abs(param_nfxp_est_array(:,:,4)-param_SLC_est_array(:,:,1)),[],'all');
param_diff_6400_nfxp_numer_EVFI=max(abs(param_nfxp_est_array(:,:,5)-param_SLC_est_array(:,:,1)),[],'all');

param_diff_1600_nfxp_analy=max(abs(param_nfxp_est_array(:,:,3)-param_SLC_est_array(:,:,3)),[],'all');

results_param_est_table=...
    [results_param_est{1};...
    results_param_est{2};...
    results_param_est{3}];

results_comp_table=...
    [results_comp{1};...
    results_comp{2};...
    results_comp{3}];

writematrix(round(results_comp_table,1),append(save_path,"dynamic_game_unobs_hetero_table_comp.csv"))

writematrix(round(results_param_est_table,3),append(save_path,"dynamic_game_unobs_hetero_table_param_est.csv"))
param_diff_table=[param_diff_6400_nfxp_analy;param_diff_6400_SLC_numer;...
    param_diff_6400_nfxp_numer;param_diff_1600_nfxp_analy];

writematrix(round([speed_test_info.comp_time_eval_Phi_func;...
    speed_test_info.comp_time_eval_obj_func],4),append(save_path,"dynamic_game_unobs_hetero_table_comp_time_obj_G.csv"))

writematrix(round(results_comp_table,1),append(save_path,"dynamic_game_unobs_hetero_table_comp.csv"))

%% EVFI comparison (Analytical derivatives)
table=round([results_comp{1};results_comp{4}],1);
table2=[0;param_diff_6400_nfxp_analy;...
    param_diff_6400_SLC_analy_EVFI;param_diff_6400_nfxp_analy_EVFI];
table=[table,round(table2,8)];

writematrix(table,append(save_path,"dynamic_game_unobs_hetero_table_comp_analy_diff_EVFI.csv"))

%% EVFI comparison (Numerical derivatives)
table=round([results_comp{2};results_comp{5}],1);
table2=[param_diff_6400_SLC_numer;param_diff_6400_nfxp_numer;...
    param_diff_6400_SLC_numer_EVFI;param_diff_6400_nfxp_numer_EVFI];
table=[table,round(table2,8)];

writematrix(table,append(save_path,"dynamic_game_unobs_hetero_table_comp_numer_diff_EVFI.csv"))
