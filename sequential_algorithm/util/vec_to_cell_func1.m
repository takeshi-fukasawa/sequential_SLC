function X_cell=vec_to_cell_func1(X_vec,cell_size)

    n_var_type=size(cell_size,2);

    X_cell=[];
    id=1;
    for i=1:n_var_type
        id_end=id+prod(cell_size{i})-1;
        X_cell{i}=reshape(X_vec(id:id_end),cell_size{i});
        id=id_end+1;
    end

end
