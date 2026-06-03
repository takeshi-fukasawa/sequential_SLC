%run("./addpath_preliminary.m")
save_path="C:/Users/fukas/Dropbox/sequential/output/";
%save_path="C:/Users/fukas/Dropbox/sequential/output/temp_0520_no_constraints/";

nrepli=20;
sim_id_max=1;

n_param=3;

%load(append(save_path,"dynamic_BLP_seq_results_SLC_only.mat"),...
%    'obj_val_SLC_all','param_SLC_all','iter_info_SLC','iter_SLC_all','times_SLC_all',...
%    'count_G_eval_SLC_all','count_obj_eval_SLC_all')

load(append(save_path,"dynamic_BLP_seq_results.mat"));%,...
        %'obj_val_nfxp_all','param_nfxp_all','iter_info_nfxp','iter_nfxp_all','times_nfxp_all',...
        %'count_G_eval_nfxp_all','count_obj_eval_nfxp_all',...
        %'obj_val_SLC_all','param_SLC_all','iter_info_SLC','iter_SLC_all','times_SLC_all',...
        %'count_G_eval_SLC_all','count_obj_eval_SLC_all','iter_info_SLC_all')

temp=load(append(save_path,"dynamic_BLP_seq_results_SLC_only_spectral.mat"));
obj_val_SLC_spectral_all=temp.obj_val_SLC_all;
param_SLC_spectral_all=temp.param_SLC_all;
iter_info_SLC_spectral_all=temp.iter_info_SLC_all;
times_SLC_spectral_all=temp.times_SLC_all;
count_G_eval_SLC_spectral_all=temp.count_G_eval_SLC_all;
count_obj_eval_SLC_spectral_all=temp.count_obj_eval_SLC_all;
iter_SLC_spectral_all=temp.iter_SLC_all;

temp=load(append(save_path,"dynamic_BLP_seq_results_SLC_only_SLC_loose_tol.mat"));
obj_val_SLC_loose_all=temp.obj_val_SLC_all;
param_SLC_loose_all=temp.param_SLC_all;
iter_info_SLC_loose_all=temp.iter_info_SLC_all;
times_SLC_loose_all=temp.times_SLC_all;
count_G_eval_SLC_loose_all=temp.count_G_eval_SLC_all;
count_obj_eval_SLC_loose_all=temp.count_obj_eval_SLC_all;
iter_SLC_loose_all=temp.iter_SLC_all;

for sim_id=1:sim_id_max


    if sim_id==1
        results_param_est=cell(sim_id_max,1);
        results_comp=cell(sim_id_max,1);

        param_nfxp_est_array=NaN(nrepli,n_param,3);
        param_SLC_est_array=NaN(nrepli,n_param,3);
        param_SLC_spectral_est_array=NaN(nrepli,n_param,3);
        
    end

    MAX_ITER_SLC=50;
    non_conv_SLC=(iter_SLC_all+1>=MAX_ITER_SLC|iter_SLC_all==1);
    non_conv_SLC_spectral=(iter_SLC_spectral_all+1>=MAX_ITER_SLC|iter_SLC_spectral_all==1);
    non_conv_SLC_loose=(iter_SLC_loose_all+1>=MAX_ITER_SLC|iter_SLC_loose_all==1);

    MAX_ITER_nfxp=100;%%%
    non_conv_nfxp=(iter_nfxp_all>=MAX_ITER_nfxp|iter_nfxp_all==1);
    obj_val_SLC_all(non_conv_SLC)=1e10;
    obj_val_SLC_spectral_all(non_conv_SLC_spectral)=1e10;
    obj_val_SLC_loose_all(non_conv_SLC_loose)=1e10;
    obj_val_nfxp_all(non_conv_nfxp)=1e10;

    [min_val_SLC, temp] = min(obj_val_SLC_all, [], 2);
    param_SLC_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_SLC_all_i = param_SLC_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_SLC_all_i), rows, temp);
        param_SLC_est(:,i) = param_SLC_all_i(idx);
    end

    [min_val_SLC_spectral, temp] = min(obj_val_SLC_spectral_all, [], 2);
    param_SLC_spectral_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_SLC_spectral_all_i = param_SLC_spectral_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_SLC_spectral_all_i), rows, temp);
        param_SLC_spectral_est(:,i) = param_SLC_spectral_all_i(idx);
    end

    [min_val_SLC_loose, temp] = min(obj_val_SLC_loose_all, [], 2);
    param_SLC_loose_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_SLC_loose_all_i = param_SLC_loose_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_SLC_loose_all_i), rows, temp);
        param_SLC_loose_est(:,i) = param_SLC_loose_all_i(idx);
    end
    
    [min_val_nfxp, temp] = min(obj_val_nfxp_all, [], 2);
    param_nfxp_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_nfxp_all_i = param_nfxp_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_nfxp_all_i), rows, temp);
        param_nfxp_est(:,i) = param_nfxp_all_i(idx);
    end

    
    times_nfxp_sum=sum(times_nfxp_all,2);
    times_SLC_sum=sum(times_SLC_all,2);
    times_SLC_spectral_sum=sum(times_SLC_spectral_all,2);
    times_SLC_loose_sum=sum(times_SLC_loose_all,2);
    %%%%
    
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
    col_SLC_spectral=[results_time_dist_func(sum(times_SLC_spectral_all,2)),...
        mean(count_obj_eval_SLC_spectral_all(:),'omitnan'),...
        mean(count_G_eval_SLC_spectral_all(:),'omitnan'),...
        mean(iter_SLC_spectral_all(:)+1),...
        mean(1-min(non_conv_SLC_spectral,[],2))*100,mean(1-non_conv_SLC_spectral(:))*100];
    col_SLC_loose=[results_time_dist_func(sum(times_SLC_loose_all,2)),...
        mean(count_obj_eval_SLC_loose_all(:),'omitnan'),...
        mean(count_G_eval_SLC_loose_all(:),'omitnan'),...
        mean(iter_SLC_loose_all(:)+1),...
        mean(1-min(non_conv_SLC_loose,[],2))*100,mean(1-non_conv_SLC_loose(:))*100];

    
    results_comp_table=[col_SLC;col_SLC_loose;col_SLC_spectral;col_nfxp];
    results_param_est=[mean(abs(param_SLC_est));std(abs(param_SLC_est));...
            mean(abs(param_SLC_spectral_est));std(abs(param_SLC_spectral_est));...
        mean(abs(param_nfxp_est));std(abs(param_nfxp_est))];

    param_nfxp_est_array(:,:,sim_id)=param_nfxp_est;
    param_SLC_est_array(:,:,sim_id)=param_SLC_est;
    param_SLC_loose_est_array(:,:,sim_id)=param_SLC_loose_est;
    param_SLC_spectral_est_array(:,:,sim_id)=param_SLC_spectral_est;

end

param_diff_nfxp_analy=max(abs(param_nfxp_est_array(:,:,1)-param_SLC_est_array(:,:,1)),[],'all');



writematrix(round(results_comp_table,1),append(save_path,"dynamic_BLP_table_comp.csv"))
writematrix(round(results_param_est,3),append(save_path,"dynamic_BLP_table_param_est.csv"))

min_val_difference_SLC=min(min_val_SLC,[],2)-min(min_val_nfxp,[],2)
min_val_difference_SLC_spectral=min(min_val_SLC_spectral,[],2)-min(min_val_nfxp,[],2)
