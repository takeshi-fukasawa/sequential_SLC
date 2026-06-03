function [output,other_vars]=v_psteady_update_func(v,psteady_1,param,...
middle_vars,numz_exo_vec,bdisc,aval,mstate)
  
global VFI_spec

  pi_param0=param(end);
  nplayer = size(v, 2);
  n_state=size(v,1);
  numa = 2^nplayer;
  numz_exo=n_state/numa;
  nalt=size(v,3);

  theta=param(1:end-1);


   %if isempty(middle_vars)==1 % If middle vars are empty, compute using functions 

        middle_vars=...
            compute_e_u_from_v_pi_func_unobs_hetero(v,psteady_1, pi_param0,numz_exo_vec, bdisc,...
        aval,mstate);

   %end
   e0=middle_vars{1};
   e1=middle_vars{2};
   u0=middle_vars{3};
   u1=middle_vars{4};
   ptrana=middle_vars{5};

    % ----------------------------
    %  (e) Updating probabilities
    % ----------------------------
    v_updated=v;
    
    for i = 1:nplayer
      v_updated(:,i,2) = u1((i-1)*n_state+1:i*n_state,:)*theta + e1((i-1)*n_state+1:i*n_state,:);
      v_updated(:,i,1) = u0((i-1)*n_state+1:i*n_state,:)*theta + e0((i-1)*n_state+1:i*n_state,:);      
    end
    

    %%% ptrana:n_state*numa;
    %%% ptrana(i,j)=Prob(State j|State i)

    ptranz_exo=compute_ptranz_exo_func(pi_param0,numz_exo_vec);
    ptranz_exo_kron=(kron(ptranz_exo, ones(numa,numa)));
  
    ptran = ptranz_exo_kron .* (kron(ones(1,numz_exo), ptrana));%n_state*n_state

    
    %psteady_updated = ptran' * psteady;
    psteady_updated = ptran' * [psteady_1;1-sum(psteady_1)];
    psteady_1_updated=psteady_updated(1:end-1);

    other_vars=[];

    if VFI_spec==1 % Relative VFI by Bray (2019)
        %other_vars.v_updated_original=v_updated;
        v_updated=v_updated-v_updated(1,:,1);% If different values are subtracted across alternatives, strange outcome...
    elseif VFI_spec==2 % Endogenous VFI by Bray (2019)
        %other_vars.v_updated_original=v_updated;
        v_updated_reshaped=reshape(v_updated,numa,numz_exo,nplayer,nalt);
        v_updated_reshaped=v_updated_reshaped-mean(v_updated_reshaped(:,:,:,1),1);% If different values are subtracted across alternatives, strange outcome...
        v_updated=reshape(v_updated_reshaped,n_state,nplayer,nalt);
    end
    

    output={v_updated,psteady_1_updated};

end
