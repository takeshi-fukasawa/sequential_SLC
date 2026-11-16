function [dG_dparam]=...
    compute_dG_dparam_func_given_v_psteady(v,psteady_1,param,...
    middle_vars,numz_exo_vec,bdisc,aval,mstate)

    psteady=[psteady_1;1-sum(psteady_1)];
    [n_state,nplayer,nalt]=size(v);
    numa=size(aval,1);
    n_param=size(param(:),1);
    naction=nalt;

    eulerc = 0.5772;
    myzero = 1e-16;
    numz_exo=prod(numz_exo_vec(:));
    kparam=size(mstate,2)+2;
    
    if size(numz_exo_vec,2)==1
        theta_rn=param(nplayer+2);
    else
        theta_rn=param(nplayer+3);
    end

    pi_param0=param(end);
    theta_initial=param(1:end-1);

    [ptranz_exo,grad_ptranz_exo]=compute_ptranz_exo_func(pi_param0,numz_exo_vec);
    ptranz_exo_kron=(kron(ptranz_exo, ones(numa,numa)));
    grad_ptranz_exo_kron=(kron(grad_ptranz_exo, ones(numa,numa)));
 
    % ----------------------------------------------------------------
    % Initialize storage for flow utilities, continuation values, etc.
    % ----------------------------------------------------------------
    u = zeros(n_state, kparam, nplayer, naction);
    fv = zeros(n_state, nplayer, naction);
    pchoice = zeros(n_state, nplayer);
    ssurplus = zeros(n_state, nplayer);
        

    % ------------------------------------------------
    % Matrices of choice probabilities,
    % social surplus, and
    % derivatives of choice probabilities.
    % ------------------------------------------------
    for i = 1:nplayer
        vdiff = v(:,i,2) - v(:,i,1);
        expvdiff = exp(vdiff);
        pchoice(:,i) = expvdiff ./ (1 + expvdiff);
        ssurplus(:,i) = log(expvdiff + 1) + v(:,i,1) + eulerc;
    end

    % -----------------------------------------------------------
    % Matrix of transition probabilities Pr(a[t]|s[t],a[t-1])
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

    % Reset Jacobian index counter
    idx = 0;
    for i = 1:nplayer
        % --------------------------------------------
        % Matrices Pr(a[t] | s[t], a[t-1], ai[t])
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

        % --------------------------------------------
        % Transition probability matrix Pr(s[t+1] | s[t], a[t-1], ai[t])
        % --------------------------------------------
        Fi0 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai0);
        Fi1 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai1);

        Fi0_grad = kron(grad_ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai0);
        Fi1_grad = kron(grad_ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai1);

        % ------------------------------------
        % Computing hi = E(ln(Sum(aj)+1))
        % ------------------------------------
        hi = aval;
        hi(:,i) = ones(numa, 1);
        hi = ptranai1 * log(sum(hi, 2));

        % ---------------------------
        % Creating U and FV
        % ---------------------------

        % Player i chooses 0
        u(:,:,i,1) = zeros(n_state, kparam); % u is initialized with zeros
        fv(:,i,1) = bdisc * Fi0 * ssurplus(:,i);      % Player i chooses 0
        fv_grad(:,i,1) = bdisc * Fi0_grad * ssurplus(:,i);    

        % Player i chooses 1
        u1i = zeros(1, nplayer);
        u1i(1,i) = 1; % Dummy for FC for player i

        if size(mstate,2)==nplayer+1
            u(:,:,i,2) = [ repmat(u1i, n_state, 1), mstate(:,1), (-hi), (mstate(:,i+1)-1) ];
        else
            u(:,:,i,2) = [ repmat(u1i, n_state, 1), mstate(:,1),mstate(:,nplayer+2:end), (-hi), (mstate(:,i+1)-1) ];
        end
        fv(:,i,2) = bdisc * Fi1 * ssurplus(:,i);      % Player i chooses 1
        fv_grad(:,i,2) = bdisc * Fi1_grad * ssurplus(:,i);    
    end

    % -------------------
    % Reshape utility function (u), future values (fv), and value function (v) to conform
    % -------------------
    fv_grad_vec=fv_grad(:);% Derivative wrt pi(state transition param)
    uperm = permute(u, [ 1, 3, 4, 2 ]);
    uvec = reshape(uperm, [n_state * nplayer * nalt, kparam]);


    %% psteady
    grad_ptran = grad_ptranz_exo_kron .* (kron(ones(1,numz_exo), ptrana));%n_state*n_state

    grad_psteady_updated = grad_ptran' * psteady;
    grad_psteady_1_updated=grad_psteady_updated(1:end-1);

    %%% G(Y;param)=Y-Phi(Y,param)

    term_v = [ -uvec,-fv_grad_vec ];%n_v*n_param
    term_p = [ zeros(n_state-1,n_param-1),...
            -grad_psteady_1_updated];%n_psteady_1*n_param
    dG_dparam=[term_v;term_p];
    
end

