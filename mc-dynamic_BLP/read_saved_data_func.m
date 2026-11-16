function [obj_val_SLC_all,param_SLC_all,iter_info_SLC_all,times_SLC_all,...
    count_G_eval_SLC_all,count_obj_eval_SLC_all,iter_SLC_all,feval_SLC_all,div_SLC_all]=...
    read_saved_data_func(filename,save_path)

    temp=load(append(save_path,filename));

    obj_val_SLC_all=temp.obj_val_SLC_all;
    param_SLC_all=temp.param_SLC_all;
    iter_info_SLC_all=temp.iter_info_SLC_all;
    times_SLC_all=temp.times_SLC_all;
    count_G_eval_SLC_all=temp.count_G_eval_SLC_all;
    count_obj_eval_SLC_all=temp.count_obj_eval_SLC_all;
    iter_SLC_all=temp.iter_SLC_all;

    [n_sim,nstart]=size(iter_SLC_all);
    for i=1:n_sim
        for id=1:nstart
            iter_info_SLC=iter_info_SLC_all{i,id};
            %iter_SLC_all(i,id)=iter_info_SLC.n_iter;
            feval_SLC_all(i,id)=iter_info_SLC.feval;
            div_SLC_all(i,id)=isnan(iter_info_SLC.count_obj_eval);
        end
    end
end
