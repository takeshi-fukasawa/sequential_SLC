function [param_SLC_est,min_val_SLC]=obtain_param_SLC_func(param_SLC_all,obj_val_SLC_all)
    [nrepli,nstart,n_param]=size(param_SLC_all);
    [min_val_SLC, temp] = min(obj_val_SLC_all, [], 2);
    param_SLC_est = NaN(nrepli, n_param);
    rows = (1:nrepli)'; 
    for i = 1:n_param
        param_SLC_all_i = param_SLC_all(:,:,i); % nrepli * nstart_seq
        idx = sub2ind(size(param_SLC_all_i), rows, temp);
        param_SLC_est(:,i) = param_SLC_all_i(idx);
    end
end
