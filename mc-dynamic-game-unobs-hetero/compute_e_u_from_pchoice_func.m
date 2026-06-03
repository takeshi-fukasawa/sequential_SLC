function out=compute_e_u_from_pchoice_func...
    (pchoice, ptranz_exo, bdisc,...
    aval,mstate,ptranz_exo_kron)

  % ---------------
  % Some constants
  % ---------------
  eulerc = 0.5772;
  myzero = 1e-16;
  nplayer = size(aval, 2);
  numa = 2^nplayer;
  n_state=size(pchoice,1);
  numz_exo=n_state/numa;
  kparam=size(mstate,2)+2;

    u0 = zeros(n_state*nplayer,kparam);
    u1 = zeros(n_state*nplayer,kparam);
    e0 = zeros(n_state*nplayer,1);
    e1 = zeros(n_state*nplayer,1);

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

    % -----------------------
    %  (b) Storing I-b*F
    % -----------------------
    i_bf = kron(ptranz_exo, ones(numa, numa)).*kron(ones(1,numz_exo), ptrana);
    i_bf = eye(n_state) - bdisc*i_bf;

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


      % ------------------------------
      %  (c.4) Creating sumu and sume
      % ------------------------------
      ppi1 = kron(ppi, ones(1, kparam));
      ppi0 = 1 - ppi1;
      sumu = ppi0.*umat0 + ppi1.*umat1;
      sume = (1-ppi).*(eulerc-log(1-ppi)) + ppi.*(eulerc-log(ppi));

      % -------------------
      %  (c.5) Creating ww
      % -------------------
      rhs = [ sumu, sume ];
      ww = i_bf \ rhs;

      % ----------------------------------
      %  (c.6) Creating utilda and etilda
      % ----------------------------------
      %%ptranai0 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai0);

      ptranai0=reshape(reshape(ptranz_exo_kron,n_state*numa,numz_exo).*...
          reshape(ptranai0,n_state*numa,1),n_state,numa*numz_exo);

      utilda0 = umat0 + bdisc*(ptranai0*ww(:,1:kparam));
      etilda0 = bdisc*(ptranai0*ww(:,kparam+1));

      %%ptranai1 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai1);
      
      ptranai1=reshape(reshape(ptranz_exo_kron,n_state*numa,numz_exo).*...
          reshape(ptranai1,n_state*numa,1),n_state,numa*numz_exo);

      utilda1 = umat1 + bdisc*(ptranai1*ww(:,1:kparam));
      etilda1 = bdisc*(ptranai1*ww(:,kparam+1));

      % -------------------------------------------
      %  (c.7) Creating observations uobs and eobs
      % -------------------------------------------
      u0((i-1)*n_state+1:i*n_state,:) = utilda0;
      u1((i-1)*n_state+1:i*n_state,:) = utilda1;
      e0((i-1)*n_state+1:i*n_state,:) = etilda0;
      e1((i-1)*n_state+1:i*n_state,:) = etilda1;

    end

    out={e0,e1,u0,u1};

end
