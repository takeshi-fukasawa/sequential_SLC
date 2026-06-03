function out=compute_e_u_from_v_func...
    (v, ptranz_exo, bdisc,...
    aval,mstate,ptranz_exo_kron)

  % ---------------
  % Some constants
  % ---------------
  eulerc = 0.5772;
  myzero = 1e-16;
  nplayer = size(aval, 2);
  numa = 2^nplayer;
  n_state=size(v,1);
  numz_exo=n_state/numa;
  kparam=size(mstate,2)+2;

    u0 = zeros(n_state*nplayer,kparam);
    u1 = zeros(n_state*nplayer,kparam);
    e0 = zeros(n_state*nplayer,1);
    e1 = zeros(n_state*nplayer,1);

    pchoice=zeros(n_state, nplayer);
    for i = 1:nplayer
        vdiff = v(:,i,2) - v(:,i,1);
        expvdiff = exp(vdiff);
        pchoice(:,i) = expvdiff ./ (1 + expvdiff);
    end
    V=log(sum(exp(v),3))+eulerc;

        % -----------------------------------------------------------
    % (a) Matrix of transition probabilities Pr(a[t]|s[t],a[t-1])
    % -----------------------------------------------------------
    ptrana = ones(n_state, numa);
    for i = 1:nplayer;
      mi = aval(:,i)';
      ppi = pchoice(:,i);
      ppi1 = repmat(ppi, 1, numa);
      ppi0 = 1 - ppi1;
      mi1 = repmat(mi, n_state, 1);
      mi0 = 1 - mi1;
      ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
    end

    % -----------------------------------------
    %  (c) Construction of explanatory variables
    % -----------------------------------------
    
    for i = 1:nplayer
      % --------------------------------------------
      %  (c.1) Matrices Pr(a[t] | s[t],a[t-1], ai[t])
      % --------------------------------------------
      mi = aval(:,i)';
      ppi = pchoice(:,i);
      ppi= (ppi>=myzero).*(ppi<=(1-myzero)).*ppi ...
         + (ppi<myzero).*myzero ...
         + (ppi>(1-myzero)).*(1-myzero);
      ppi1 = repmat(ppi, 1, numa);
      ppi0 = 1 - ppi1;
      mi1 = repmat(mi, n_state, 1);
      mi0 = 1 - mi1;
      ptrani = ((ppi1 .^ mi1) .* (ppi0 .^ mi0));
      ptranai0 = ptrana .* (mi0 ./ ptrani);
      ptranai1 = ptrana .* (mi1 ./ ptrani);

      % ------------------------------------
      %  (c.2) Computing hi = E(ln(Sum(aj)+1))
      % ------------------------------------
      hi = aval;
      hi(:,i) = ones(numa, 1);
      hi = ptranai1 * log(sum(hi, 2));

      % ---------------------------
      %  (c.3) Creating U0 and U1
      % ---------------------------
      umat0 = zeros(n_state, kparam);
      umat1 = eye(nplayer);
      umat1 = umat1(i,:);
      if size(mstate,2)==nplayer+1
        umat1 = [ repmat(umat1, n_state, 1), mstate(:,1), (-hi), (mstate(:,i+1)-1) ];
      else
        umat1 = [ repmat(umat1, n_state, 1), mstate(:,1),mstate(:,nplayer+2:end), (-hi), (mstate(:,i+1)-1) ];
      end

      % ----------------------------------
      %  (c.6) Creating utilda and etilda
      % ----------------------------------
      utilda0=umat0;
      utilda1=umat1;

      %I_numa=ones(numa,numa);
      %%ptranz_exo_kron=(kron(ptranz_exo, I_numa));


      %ptranz_exo_kron=sparse(ptranz_exo_kron);
      %ptranai0=sparse(ptranai0);
      %ptranai1=sparse(ptranai1);
   
      %ones_numz_exo=ones(1,numz_exo);

      %%%%ptranai0_temp =  ptranz_exo_kron.* kron(ones_numz_exo, ptranai0);

      ptranai0_temp=reshape(reshape(ptranz_exo_kron,n_state*numa,numz_exo).*...
          reshape(ptranai0,n_state*numa,1),n_state,numa*numz_exo);

      etilda0=bdisc*(ptranai0_temp*V(:,i));

      ptranai1_temp=reshape(reshape(ptranz_exo_kron,n_state*numa,numz_exo).*...
          reshape(ptranai1,n_state*numa,1),n_state,numa*numz_exo);

      %%%%ptranai1_temp = ptranz_exo_kron .* kron(ones_numz_exo, ptranai1);
      etilda1=bdisc*(ptranai1_temp*V(:,i));
   
      % -------------------------------------------
      %  (c.7) Creating observations uobs and eobs
      % -------------------------------------------
      u0((i-1)*n_state+1:i*n_state,:) = utilda0;
      u1((i-1)*n_state+1:i*n_state,:) = utilda1;
      e0((i-1)*n_state+1:i*n_state,:) = etilda0;
      e1((i-1)*n_state+1:i*n_state,:) = etilda1;

          
    end

    out={e0,e1,u0,u1,ptrana};

end
