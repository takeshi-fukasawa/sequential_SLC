function [x_updated_cell,Delta_x_cell,G_vec,G_cell,other_vars]=...
    compute_z_func_given_v_psteady_unobs_hetero(v,psteady_1,param,...
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
    pchoice_est = zeros(n_state, nplayer);
    ssurplus = zeros(n_state, nplayer);
    pchoice_dv = zeros(n_state, naction, naction, nplayer); % Derivative of Lambda(x,a,i) wrt v(x,a',i)
    
    % --------------------------
    % Initialize Jacobian matrix
    % --------------------------
    
    % Jacobian dimension: sparse (jacdim x jacdim)
    jacdim = n_state * naction * nplayer;
    
    % Number of nonzeros in Jacobian
    % There are (naction * nplayer) x (naction * nplayer) diagonal blocks of dimension (n_state x n_state).
    jacnz = n_state*(naction*nplayer)^2 - 4*n_state*nplayer + nplayer*4*n_state^2;
    
    % Triplets for sparse jacobian
    I = zeros(jacnz, 1);
    J = zeros(jacnz, 1);
    X = zeros(jacnz, 1);


    %%% compute the matrices of choice probabilities,
    for i = 1:nplayer
        vdiff = v(:,i,2) - v(:,i,1);
        expvdiff = exp(vdiff);
        pchoice(:,i) = expvdiff ./ (1 + expvdiff);
    end

    % ------------------------------------------------
    % Matrices of choice probabilities,
    % social surplus, and
    % derivatives of choice probabilities.
    % ------------------------------------------------
    for i = 1:nplayer
        vdiff = v(:,i,2) - v(:,i,1);
        expvdiff = exp(vdiff);
        %pchoice(:,i) = expvdiff ./ (1 + expvdiff);
        ssurplus(:,i) = log(expvdiff + 1) + v(:,i,1) + eulerc;

        % Derivative of Lambda_^i(:,a^i=2) with respect to v^i(:,a^i=2)
        pchoice_dv(:,2,2,i) = pchoice(:,i) .* (1 - pchoice(:,i));
        % Derivative of Lambda_^i(:,a^i=2) with respect to v^i(:,a^i=1)
        pchoice_dv(:,2,1,i) = -pchoice_dv(:,2,2,i);
        % Derivative of Lambda_^i(:,a^i=1) with respect to v^i(:,a^i=1)
        pchoice_dv(:,1,1,i) = pchoice_dv(:,2,2,i);
        % Derivative of Lambda_^i(:,a^i=1) with respect to v^i(:,a^i=2)
        pchoice_dv(:,1,2,i) = -pchoice_dv(:,1,1,i);
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

    dGp_dv_temp=zeros(n_state-1,n_state,nplayer,nalt);

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

        %%%%%%%%%%%%%%
        %% Value for computing dGv_dv
        ppi1_dvi0 = repmat(pchoice_dv(:,2,1,i), 1, numa);
        ppi1_dvi1 = repmat(pchoice_dv(:,2,2,i), 1, numa);
        ppi0_dvi0 = repmat(pchoice_dv(:,1,1,i), 1, numa);
        ppi0_dvi1 = repmat(pchoice_dv(:,1,2,i), 1, numa);

        ptrani_dvi0 = ((ppi1_dvi0 .^ mi1) .* (ppi0_dvi0 .^ mi0));
        ptrani_dvi1 = ((ppi1_dvi1 .^ mi1) .* (ppi0_dvi1 .^ mi0));

        ptrana_dvi0= ptrana ./ ptrani .*ptrani_dvi0;
        F_dvi0 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1, numz_exo), ptrana_dvi0);%n_state*n_state
        ptrana_dvi1= ptrana ./ ptrani .*ptrani_dvi1;
        F_dvi1 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1, numz_exo), ptrana_dvi1);%n_state*n_state

        dGp_dvi0_temp=-(F_dvi0').*psteady';%n_state*n_state
        dGp_dvi1_temp=-(F_dvi1').*psteady';%n_state*n_state

        %%% Omit the last row, because psteady(end) is not independent 
        dGp_dv_temp(:,:,i,1)=dGp_dvi0_temp(1:end-1,:);%(n_state-1)*n_state
        dGp_dv_temp(:,:,i,2)=dGp_dvi1_temp(1:end-1,:);%(n_state-1)*n_state

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
    
        % ------------------------------
        % Calculate the Jacobian
        % ------------------------------
        %
        % The Jacobian is a sparse matrix with the following block structure.
        % Each (ai,aj) cell represents a matrix of dimension (n_state x n_state).
        %
        % [  (a1=0,a1=0)  (a1=0,a2=0)  ...  (a1=0,aN=0) |  (a1=0,a1=1)  (a1=0,a2=1)  ...  (a1=0,aN=1) ]
        % [  (a2=0,a1=0)  (a2=0,a2=0)  ...  (a2=0,aN=0) |  (a2=0,a1=1)  (a2=0,a2=1)  ...  (a2=0,aN=1) ]
        % [      ...          ...      ...      ...     |      ...          ...      ...      ...     ]
        % [  (aN=0,a1=0)  (aN=0,a2=0)  ...  (aN=0,aN=0) |  (aN=0,a1=1)  (aN=0,a2=1)  ...  (aN=0,aN=1) ]
        % [---------------------------------------------|---------------------------------------------]
        % [  (a1=1,a1=0)  (a1=1,a2=0)  ...  (a1=1,aN=0) |  (a1=1,a1=1)  (a1=1,a2=1)  ...  (a1=1,aN=1) ]
        % [  (a2=1,a1=0)  (a2=1,a2=0)  ...  (a2=1,aN=0) |  (a2=1,a1=1)  (a2=1,a2=1)  ...  (a2=1,aN=1) ]
        % [      ...          ...      ...      ...     |      ...          ...      ...      ...     ]
        % [  (aN=1,a1=0)  (aN=1,a2=0)  ...  (aN=1,aN=0) |  (aN=1,a1=1)  (aN=1,a2=1)  ...  (aN=1,aN=1) ]

        if size(numz_exo_vec,2)==1
            theta_rn=param(nplayer+2);
        else
            theta_rn=param(nplayer+3);
        end

        for j = 1:nplayer
            if j == i
                % ------------------------------
                % Diagonal Jacobian blocks
                % ------------------------------

                % First, store diagonal elements of these blocks as vectors (for clarity)
                jac_i0_i0 = eye(n_state,n_state) - bdisc * Fi0 .* repmat((1-ppi)', n_state, 1);
                jac_i1_i1 = eye(n_state,n_state) - bdisc * Fi1 .* repmat(ppi', n_state, 1);
                jac_i0_i1 = - bdisc * Fi0 .* repmat(ppi', n_state, 1);
                jac_i1_i0 = - bdisc * Fi1 .* repmat((1-ppi)', n_state, 1);

                % Calculate row and column indices
                jac_i0_row = repmat([(i-1)*n_state+1:(i-1)*n_state+n_state]', 1, n_state);
                jac_i0_col = repmat([(i-1)*n_state+1:(i-1)*n_state+n_state], n_state, 1);
                jac_i1_row = jac_i0_row + jacdim/2;
                jac_i1_col = jac_i0_col + jacdim/2;

                % Store values for sparse Jacobian
                I(idx+1:idx+4*n_state^2) = [ jac_i0_row(:); jac_i0_row(:); jac_i1_row(:); jac_i1_row(:); ];
                J(idx+1:idx+4*n_state^2) = [ jac_i0_col(:); jac_i1_col(:); jac_i0_col(:); jac_i1_col(:); ];
                X(idx+1:idx+4*n_state^2) = [ jac_i0_i0(:);  jac_i0_i1(:);  jac_i1_i0(:);  jac_i1_i1(:);  ];
                idx = idx + 4*n_state^2;
            else%j not equal to i
                % --------------------------------------------
                % Matrices Pr(a[t] | s[t], a[t-1], ai[t], aj[t])
                % --------------------------------------------
                mj = aval(:,j)';
                ppj = pchoice(:,j);
                ppj= (ppj>=myzero).*(ppj<=(1-myzero)).*ppj ...
                     + (ppj<myzero).*myzero ...
                     + (ppj>(1-myzero)).*(1-myzero);
                ppj1 = repmat(ppj, 1, numa);
                ppj0 = 1 - ppj1;
                mj1 = repmat(mj, n_state, 1);
                mj0 = 1 - mj1;
                ptranj = ((ppj1 .^ mj1) .* (ppj0 .^ mj0));

                % --------------------------------------------
                % Derivatives of transition probability Pr(s[t+1] | s[t], a[t-1], ai[t])
                % --------------------------------------------
                ppj1_dvj0 = repmat(pchoice_dv(:,2,1,j), 1, numa);
                ppj1_dvj1 = repmat(pchoice_dv(:,2,2,j), 1, numa);
                ppj0_dvj0 = repmat(pchoice_dv(:,1,1,j), 1, numa);
                ppj0_dvj1 = repmat(pchoice_dv(:,1,2,j), 1, numa);

                ptranj_dvj0 = ((ppj1_dvj0 .^ mj1) .* (ppj0_dvj0 .^ mj0));
                ptranj_dvj1 = ((ppj1_dvj1 .^ mj1) .* (ppj0_dvj1 .^ mj0));

                ptranai0_dvj0 = ptranai0 ./ ptranj .* ptranj_dvj0;
                ptranai0_dvj1 = ptranai0 ./ ptranj .* ptranj_dvj1;
                ptranai1_dvj0 = ptranai1 ./ ptranj .* ptranj_dvj0;
                ptranai1_dvj1 = ptranai1 ./ ptranj .* ptranj_dvj1;

                Fi0_dvj0 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai0_dvj0);
                Fi0_dvj1 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai0_dvj1);
                Fi1_dvj0 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai1_dvj0);
                Fi1_dvj1 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai1_dvj1);


                % ------------------------------
                % Off-diagonal Jacobian blocks
                % ------------------------------

                % Choice ai = 0 (flow utility is zero, future value term only)
                jac_i0_j0_val = -bdisc*Fi0_dvj0*ssurplus(:,i);
                jac_i0_j1_val = -bdisc*Fi0_dvj1*ssurplus(:,i);

                % Used for calculating derivative of expected log number of rivals given ai=1
                hitmp = aval;
                hitmp(:,i) = ones(numa, 1);

                utmp = -ptranai1_dvj0 * log(sum(hitmp, 2)); % Appears in u with a negative sign

                jac_i1_j0_val = -utmp.*theta_rn - bdisc*Fi1_dvj0*ssurplus(:,i);


                utmp = -ptranai1_dvj1 * log(sum(hitmp, 2)); % Appears in u with a negative sign
                jac_i1_j1_val = -utmp.*theta_rn - bdisc*Fi1_dvj1*ssurplus(:,i);

                % Calculate row and column offsets
                jac_i0_row = [(i-1)*n_state+1:(i-1)*n_state+n_state]';
                jac_i1_row = jac_i0_row + jacdim/2;
                jac_j0_col = [(j-1)*n_state+1:(j-1)*n_state+n_state]';
                jac_j1_col = jac_j0_col + jacdim/2;

                % Store values for sparse Jacobian
                I(idx+1:idx+4*n_state) = [ jac_i0_row;    jac_i0_row;    jac_i1_row;    jac_i1_row;    ];
                J(idx+1:idx+4*n_state) = [ jac_j0_col;    jac_j1_col;    jac_j0_col;    jac_j1_col;    ];
                X(idx+1:idx+4*n_state) = [ jac_i0_j0_val; jac_i0_j1_val; jac_i1_j0_val; jac_i1_j1_val; ];
                idx = idx + 4*n_state;
            end            
        end%j = 1:nplayer
    end%i = 1:nplayer

    jac_Gvv = sparse(I, J, X, jacdim, jacdim, jacnz);
    jac_Gvv=-full(jac_Gvv);%n_v*n_v%%%%%%


    % -------------------
    % Reshape utility function (u), future values (fv), and value function (v) to conform
    % -------------------
    vvec = v(:);
    fvvec = fv(:);
    fv_grad_vec=fv_grad(:);% Derivative wrt pi(state transition param)
    uperm = permute(u, [ 1, 3, 4, 2 ]);
    uvec = reshape(uperm, [n_state * nplayer * nalt, kparam]);


    %% psteady
    ptran = ptranz_exo_kron .* (kron(ones(1,numz_exo), ptrana));%n_state*n_state

    grad_ptran = grad_ptranz_exo_kron .* (kron(ones(1,numz_exo), ptrana));%n_state*n_state

    psteady_updated = ptran' * psteady;
    psteady_1_updated=psteady_updated(1:end-1);

    grad_psteady_updated = grad_ptran' * psteady;
    grad_psteady_1_updated=grad_psteady_updated(1:end-1);

    A=ptran';

    jac_Gpp=eye(n_state-1,n_state-1)-A(1:end-1,1:end-1)+A(1:end-1,end);%n_state*n_state
    n_v=size(v(:),1);%n_v=n_state*player*2
    jac_Gvp=zeros(n_v,n_state-1);% dGv_dp

    jac_Gpv=reshape(dGp_dv_temp,n_state-1,n_v);

    jac=sparse([jac_Gvv,jac_Gvp;jac_Gpv,jac_Gpp]);

    %% Solve linear equation
    G_v=vvec-(uvec*theta_initial+fvvec);%n_v*1
    G_p=psteady_1-psteady_1_updated;%n_psteady_1*1
    G_vec=[G_v(:);G_p(:)];
    G_cell={reshape(G_v,size(v));reshape(G_p,size(psteady_1))};

    rhs_Gv = [ G_v ];%n_v*(n_param+1)
    rhs_Gp = [-G_p];%n_psteady_1*(n_param+1)

    rhs=[rhs_Gv;rhs_Gp];

    ww = jac \ rhs;

    z_vec=ww(:,end);  
    z_cell{1}=reshape(z_vec(1:n_v,:),size(v));
    z_cell{2}=reshape(z_vec(n_v+1:end,:),size(psteady_1));

    x_updated_cell{1}=v+reshape(z_cell{1},size(v));;
    x_updated_cell{2}=psteady_1+reshape(z_cell{2},size(psteady_1));

    Delta_x_cell{1}=x_updated_cell{1}-v;
    Delta_x_cell{2}=x_updated_cell{2}-psteady_1;

    other_vars=[];
end

