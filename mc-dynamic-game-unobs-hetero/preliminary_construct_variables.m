% ---------------
% Some constants
% ---------------
eulerc = 0.5772;
myzero = 1e-16;
nobs = size(aobs, 1);
nplayer = size(aobs, 2);
naction = 2;
numa = 2^nplayer;
numz_exo = size(z_exoval, 1);%%%
n_state = numz_exo*numa;

% ----------------------------------------------
% Matrix with values of states of (s[t],a[t-1])
% ----------------------------------------------
aval = zeros(numa, nplayer);
for i = 1:nplayer
    aval(:,i) = kron(ones(2^(i-1),1), kron([ 0; 1 ], ones(2^(nplayer-i),1)));
end

mstate=compute_mstate_func(z_exoval,n_state,nplayer,numa,numz_exo,aval);
kparam = size(mstate,2)+2;

critconv = (1e-2)*(1/kparam);
times = [];
varb = zeros(kparam, kparam);
llike = -1e300;
% ------------------------------------------------
% Matrix with observed indexes of state variables
% ------------------------------------------------

[indobs_orig,indobs_a_1] = indobs_func(aobs_1,zexo_obs,z_exoval);
indobs=indobs_orig;
for i = 2:nplayer
    indobs = [ indobs; indobs_orig + (i-1)*n_state ];
end

dummy_mat=(reshape(aobs,nobs,1,nplayer)==...
    reshape(0:naction-1,1,naction,1));%nobs*naction*nplayer

dummy_count=zeros(n_state*nplayer,naction);

for i=1:(n_state*nplayer)
    dummy_count(i,:)=sum(reshape(permute(dummy_mat,[1,3,2]),[],naction).*(indobs==i),1);
end

ptranz_exo_kron=(kron(ptranz_exo, ones(numa,numa)));
