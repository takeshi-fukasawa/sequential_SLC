x_cell=x_initial_cell;
x_cell{1}=x_initial_cell{1}*1.1+rand(1);
x_cell{2}=x_initial_cell{2}*0.9+rand(1);

%%%%%%%
    %param=param_k;
    eps_val=1e-4;
    id=size(param_initial(:),1);
    %id=2;
    [obj_val,dobj_dparam,dobj_dX_cell]=obj_func(param_initial,x_cell{:},other_input_obj_func_test{:});
    param_initial_eps=param_initial;
    param_initial_eps(id)=param_initial_eps(id)+eps_val;
    
    obj_val_eps=obj_func(param_initial_eps,x_cell{:},other_input_obj_func_test{:});

    grad_numer_param=(obj_val_eps-obj_val)./eps_val
    grad_analy_param=dobj_dparam(id)

    %%%%%%%
    var_id=2;
    id=size(dobj_dX_cell{var_id}(:),1);
    
    id=8;
    x_cell_eps=x_cell;
    x_cell_eps{var_id}(id)=x_cell_eps{var_id}(id)+eps_val;
    obj_val_eps=obj_func(param_initial,x_cell_eps{:},...
        other_input_obj_func_test{:});

    grad_numer=(obj_val_eps-obj_val)./eps_val
    grad_analy=dobj_dX_cell{var_id}(id)

    temp=0;