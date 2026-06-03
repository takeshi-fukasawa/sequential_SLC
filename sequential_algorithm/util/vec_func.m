function out_vec=vec_func(X_vec,cell_func,cell_size,varargin)

    X_cell=vec_to_cell_func1(X_vec,cell_size);
    out_cell=cell_func(X_cell{:},varargin{:});

    out_vec=cell_to_vec_func(out_cell);

end
