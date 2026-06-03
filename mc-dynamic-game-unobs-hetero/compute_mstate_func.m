function mstate=compute_mstate_func(z_exoval,n_state,nplayer,numa,numz_exo,aval)
    mstate = zeros(n_state, nplayer+size(z_exoval,2));
   temp=kron(z_exoval, ones(numa, 1));
    mstate(:,1) = temp(:,1);
    mstate(:,2:nplayer+1) = kron(ones(numz_exo, 1), aval);

    if size(z_exoval,2)>=2
       mstate(:,nplayer+2:end) = temp(:,2:end);
    end
end