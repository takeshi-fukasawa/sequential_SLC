function [obj_val,dobj_dparam,dobj_dX_cell,likelihood_vec]=obj_func_given_v_psteady_unobs_hetero(...
    param,v,psteady,...
    middle_vars, numz_exo_vec,bdisc,...
        aval,mstate,dummy_aobs,dummy_xobs_obs,indobs_state_obs)


if nargout==1
  likelihood_vec=compute_likelihood_given_v_psteady_unobs_hetero(...
    param,v,psteady,...
    middle_vars, numz_exo_vec,bdisc,...
        aval,mstate,dummy_aobs,dummy_xobs_obs,indobs_state_obs);

else% Compute gradient
    [likelihood_vec,dobj_dparam,dobj_dv,dobj_dpsteady_1]=...
    compute_likelihood_given_v_psteady_unobs_hetero(...
        param,v,psteady,...
        middle_vars, numz_exo_vec,bdisc,...
            aval,mstate,dummy_aobs,dummy_xobs_obs,indobs_state_obs);

    dobj_dX_cell={dobj_dv,dobj_dpsteady_1};
end


obj_val=-sum(log(likelihood_vec(:)))./size(likelihood_vec,1);

%%%%%%%%%%%
%%% comp time test
%pi_param_temp=param(end);
%numa=size(aval,1);
%ptranz_exo=compute_ptranz_exo_func(pi_param_temp,numz_exo_vec);
%n_state=size(ptranz_exo,1);
%A=eye(n_state,n_state)-ptranz_exo;
%b=ones(n_state,1);
%x_sol=A\b;
%%%%%%%%%%

%obj_val
end
