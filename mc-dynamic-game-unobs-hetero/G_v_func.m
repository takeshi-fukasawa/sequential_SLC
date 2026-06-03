function resid=G_v_func(v,theta,...
ptranz_exo,bdisc,aval,mstate,ptranz_exo_kron)

%%% v: vector; not array %%%

  nplayer = size(aval, 2);
  nstate=size(mstate,1);
  naction=2;
  v=reshape(v,nstate,nplayer,naction);

   [out,other_vars]=v_update_func(v,theta,[],...
    ptranz_exo,bdisc,aval,mstate,ptranz_exo_kron);

   v_updated=out{1};
    resid=-v_updated(:)+v(:);
end
