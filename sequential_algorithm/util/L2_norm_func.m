function F_norm=L2_norm_func(F_cell)

    n_vars=size(F_cell,2);

    F_norm=0;
    for i=1:n_vars
        F_norm=F_norm+sum(F_cell{i}(:).^2,'omitnan');
    end
    F_norm=sqrt(F_norm);
    
end


