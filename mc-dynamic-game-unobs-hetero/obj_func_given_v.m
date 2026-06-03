function [obj_val,grad]=obj_func_given_v(param,v,middle_vars,... 
 ptranz_exo,bdisc,...
        aval,mstate,dummy_count,ptranz_exo_kron)


if isempty(middle_vars)==1
    middle_vars=compute_e_u_from_v_func(v, ptranz_exo, bdisc,...
            aval,mstate,ptranz_exo_kron);
end

if nargout>=2
   [obj_val,grad]=log_likelihood_func0(param,dummy_count,...
    middle_vars{1:4});
else
   [obj_val]=log_likelihood_func0(param,dummy_count,...
    middle_vars{1:4});
end
end
