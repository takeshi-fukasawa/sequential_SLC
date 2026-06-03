function [merit_val_vec,G_cell]=compute_merit_val_func(X_cell,mu_val,...
obj_func,Phi_func,other_input_obj_func,other_input_Phi_func,...
    fixed_point_iter_spec,x_max_cell,x_min_cell)

        if isempty(other_input_obj_func)==0
            obj_val=obj_func(X_cell{:},other_input_obj_func{:});
        else
            obj_val=obj_func(X_cell{:});
        end

        if size(X_cell,2)==2
            temp={X_cell{2}};
        else
            temp=[];
            for i=1:size(X_cell,2)-1
                temp{i}=X_cell{1+i};
            end
        end
        G_cell=G_func(Phi_func,temp,fixed_point_iter_spec,X_cell{1},...
            other_input_Phi_func{:});

        G_val=0;
        for i=1:size(G_cell,2)
            G_val=G_val+sum(abs(G_cell{i}(:)),'omitnan');
        end
            
        merit_val=obj_val+mu_val*G_val;
        
        merit_val_vec=[merit_val,obj_val,G_val];
        
end