function ind=random_draw_func(prob,n_draw)
  %prob: numz_exo*1
   numz_exo=size(prob,1);

  pbuff1 = cumsum(prob);
  pbuff0 = cumsum([ 0; prob(1:numz_exo-1) ]);
  uobs = rand(n_draw, 1);
  pbuff1 = kron(pbuff1, ones(1, n_draw));
  pbuff0 = kron(pbuff0, ones(1, n_draw));
  uobs = kron(uobs, ones(1, numz_exo));
  uobs = (uobs>=(pbuff0')) .* (uobs<=(pbuff1'));%rowsum=1
  ind = uobs * [ 1:numz_exo ]';%nobs*1
end
