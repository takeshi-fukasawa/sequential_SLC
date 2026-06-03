function out=compute_e_u_from_v_pi_func_unobs_hetero...
    (v,psteady_1,pi_param0,numz_exo_vec, bdisc,...
    aval,mstate)


  nplayer = size(aval, 2);
  numa = 2^nplayer;

  ptranz_exo=compute_ptranz_exo_func(pi_param0,numz_exo_vec);
  ptranz_exo_kron=(kron(ptranz_exo, ones(numa,numa)));

  out=compute_e_u_from_v_func...
      (v, ptranz_exo, bdisc,...
      aval,mstate,ptranz_exo_kron);
end
