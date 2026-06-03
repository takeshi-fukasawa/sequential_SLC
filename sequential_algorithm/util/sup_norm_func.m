function F_norm=sup_norm_func(F_cell)

    n_vars=size(F_cell,2);

    F_norm=0;
    for i=1:n_vars
        F_norm=max(F_norm,max(abs(F_cell{i}(:))));
    end

end


