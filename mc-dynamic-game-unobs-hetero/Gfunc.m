% Gfunc.m -- Evaluates equilibrium conditions for dynamic entry/exit model
%
% FORMAT:
% [result] = Gfunc(z_exoval, ptranz_exo, theta, v, bdisc, nplayer, naction)
%
% INPUTS:
%     z_exoval    - (numz_exo x 1) vector with values of market characteristics
%     ptranz_exo  - (numz_exo x numz_exo) matrix of transition probabilities
%               of market characteristics.
%     theta   - (kparam x 1) vector of structural parameters
%     v       - (n_state*nplayer*naction x 1) value function vector
%     bdisc   - discount factor
%     nplayer - number of players
%     naction - number of actions
%
% OUTPUT:
%
%     result  - (n_state*nplayer*naction x 1) vector of equilibrium condition
%               values (i.e., the residual G(v) = v - \Gamma(v) of the fixed
%               point mapping evaluated at the input v).

function [ result ] = Gfunc(z_exoval, ptranz_exo, theta, v, bdisc, nplayer, naction)

% ---------------
% Some constants
% ---------------
eulerc = 0.5772;
numa = 2^nplayer;
numz_exo = size(z_exoval, 1);
n_state = numz_exo * numa;
kparam = size(theta,1);
sigmaeps = 1.0;


% ----------------------------------------------
% Matrix with values of states of (s[t],a[t-1])
% ----------------------------------------------
aval = zeros(numa, nplayer);
for i = 1:nplayer
    aval(:,i) = kron(ones(2^(i-1),1), kron([ 0; 1 ], ones(2^(nplayer-i),1)));
end


mstate=compute_mstate_func(z_exoval,n_state,nplayer,numa,numz_exo,aval);


% Reshape vector v into a multi-dimensional array for clarity
v = reshape(v, [n_state, nplayer, naction]);

% ------------------------------------------------
% Matrices of choice probabilities,
% social surplus, and
% derivatives of choice probabilities.
% ------------------------------------------------
pchoice = zeros(n_state, nplayer);
ssurplus = zeros(n_state, nplayer);
for i = 1:nplayer
    vdiff = v(:,i,2) - v(:,i,1);
    expvdiff = exp(vdiff);
    pchoice(:,i) = expvdiff ./ (1 + expvdiff);
    ssurplus(:,i) = log(expvdiff + 1) + v(:,i,1) + eulerc;
end

% ----------------------
% Equilibrium conditions
% ----------------------
u = zeros(n_state, kparam, nplayer, naction);
fv = zeros(n_state, nplayer, naction);
G = zeros(n_state, nplayer, naction);

% -----------------------------------------------------------
% Matrix of transition probabilities Pr(a[t]|s[t],a[t-1])
% -----------------------------------------------------------
ptrana = ones(n_state, numa);
for i = 1:nplayer
    mi = aval(:,i)';
    ppi = pchoice(:,i);
    ppi1 = repmat(ppi, 1, numa);
    ppi0 = 1 - ppi1;
    mi1 = repmat(mi, n_state, 1);
    mi0 = 1 - mi1;
    ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
end

for i = 1:nplayer
    % --------------------------------------------
    % Matrices Pr(a[t] | s[t], a[t-1], ai[t])
    % --------------------------------------------
    mi = aval(:,i)';
    ppi = pchoice(:,i);
    ppi1 = repmat(ppi, 1, numa);
    ppi0 = 1 - ppi1;
    mi1 = repmat(mi, n_state, 1);
    mi0 = 1 - mi1;
    ptrani = ((ppi1 .^ mi1) .* (ppi0 .^ mi0));
    ptranai0 = ptrana .* (mi0 ./ ptrani);
    ptranai1 = ptrana .* (mi1 ./ ptrani);
    clear mi;

    % --------------------------------------------
    % Transition probability matrix Pr(s[t+1] | s[t], a[t-1], ai[t])
    % --------------------------------------------
    Fi0 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai0);
    Fi1 = kron(ptranz_exo, ones(numa, numa)) .* kron(ones(1,numz_exo), ptranai1);

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
    %u(:,:,i,1) = zeros(n_state, kparam); % u is initialized with zeros
    fv(:,i,1) = bdisc * Fi0 * ssurplus(:,i);      % Player i chooses 0

    % Player i chooses 1
    u1i = zeros(1, nplayer);
    u1i(1,i) = 1; % Dummy for FC for player i
    if size(mstate,2)==nplayer+1
        u(:,:,i,2) = [ repmat(u1i, n_state, 1), mstate(:,1), (-hi), (mstate(:,i+1)-1) ];
    else
        u(:,:,i,2) = [ repmat(u1i, n_state, 1), mstate(:,1),mstate(:,nplayer+2:end), (-hi), (mstate(:,i+1)-1) ];
    end

    fv(:,i,2) = bdisc * Fi1 * ssurplus(:,i);      % Player i chooses 0

    % Equilibrium conditions for player i
    G(:,i,1) = v(:,i,1) - (u(:,:,i,1) * theta + bdisc * Fi0 * ssurplus(:,i));
    G(:,i,2) = v(:,i,2) - (u(:,:,i,2) * theta + bdisc * Fi1 * ssurplus(:,i));
end

% Reshape equilibrium conditions G as vector and return
result = G(:);

%%%%%%
if 1==0
[out,other_vars]=v_update_func(v,theta,...
[],ptranz_exo,bdisc,aval,mstate,ptranz_exo_kron);
v_updated=out{1};
global ratio
ratio=v_updated./v;
end

%%%%%%

% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Profit functions (a = 1):')
% profits = [ u(:,:,1,2) * theta, u(:,:,2,2) * theta, u(:,:,3,2) * theta, u(:,:,4,2) * theta, u(:,:,5,2) * theta ];
% disp(profits)
% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Choice probabilities (a = 1):')
% disp(pchoice)
% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Social surplus functions:')
% disp(ssurplus)
% fprintf('----------------------------------------------------------------------------------------\n');
% disp('Norm of equilibrium conditions:')
% disp(norm(result))
end

