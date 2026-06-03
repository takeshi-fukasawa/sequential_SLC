function [out,other_vars]=pchoice_update_func(pchoice,theta,...
middle_vars,...
ptranz_exo,bdisc,aval,mstate,ptranz_exo_kron)

  nplayer = size(pchoice, 2);
  n_state=size(pchoice,1);
  naction = 2;

   if isempty(middle_vars)==1 % If middle vars are empty, compute using functions 

    middle_vars=...
        compute_e_u_from_pchoice_func(pchoice, ptranz_exo, bdisc,...
        aval,mstate,ptranz_exo_kron);
   end
   e0=middle_vars{1};
   e1=middle_vars{2};
   u0=middle_vars{3};
   u1=middle_vars{4};

    % ----------------------------
    %  (e) Updating probabilities
    % ----------------------------
    v = zeros(n_state,nplayer,naction);
    pchoice_updated=pchoice;

    for i = 1:nplayer
      v(:,i,2) = u1((i-1)*n_state+1:i*n_state,:)*theta + e1((i-1)*n_state+1:i*n_state,:);
      v(:,i,1) = u0((i-1)*n_state+1:i*n_state,:)*theta + e0((i-1)*n_state+1:i*n_state,:);
      buff = v(:,i,2) - v(:,i,1);
      pchoice_updated(:,i) = exp(buff)./(1+exp(buff));
    end

out={pchoice_updated};
other_vars=[];

end
