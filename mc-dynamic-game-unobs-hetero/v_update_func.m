function [out,other_vars]=v_update_func(v,theta,middle_vars,...
ptranz_exo,bdisc,aval,mstate,ptranz_exo_kron)


  nplayer = size(v, 2);
  n_state=size(v,1);

   if isempty(middle_vars)==1 % If middle vars are empty, compute using functions 

    middle_vars=...
        compute_e_u_from_v_func(v, ptranz_exo, bdisc,...
        aval,mstate,ptranz_exo_kron);
   end
   e0=middle_vars{1};
   e1=middle_vars{2};
   u0=middle_vars{3};
   u1=middle_vars{4};

    % ----------------------------
    %  (e) Updating probabilities
    % ----------------------------
    v_updated=v;
    
    for i = 1:nplayer
      v_updated(:,i,2) = u1((i-1)*n_state+1:i*n_state,:)*theta + e1((i-1)*n_state+1:i*n_state,:);
      v_updated(:,i,1) = u0((i-1)*n_state+1:i*n_state,:)*theta + e0((i-1)*n_state+1:i*n_state,:);      
    end


out={v_updated};
other_vars=[];

end
