function JVP=compute_JVP_func(X_cell,d_cell,...
func,other_input_func)

    eps_val=eps^(1/3);
    for i=1:size(X_cell,2)
        X_cell_eps_plus{i}=X_cell{i}+eps_val*reshape(d_cell{i},size(X_cell{i}));
        X_cell_eps_minus{i}=X_cell{i}-eps_val*reshape(d_cell{i},size(X_cell{i}));
    end


        if isempty(other_input_func)==0
            func_val_eps_plus=func(X_cell_eps_plus{:},other_input_func{:});
            func_val_eps_minus=func(X_cell_eps_minus{:},other_input_func{:});
        else
            func_val_eps_plus=func(X_cell_eps_plus{:});
            func_val_eps_minus=func(X_cell_eps_minus{:});
        end

        if iscell(func_val_eps_plus)==0
            JVP=(func_val_eps_plus-func_val_eps_minus)./(2*eps_val);
        else
            for i=1:size(func_val_eps_plus,2)
                 JVP{i}=(func_val_eps_plus{i}-func_val_eps_minus{i})./(2*eps_val);
            end
        end


end