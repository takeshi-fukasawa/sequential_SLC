function [ptrans,grad]=compute_ptrans_func(pi_param,nums)

%% Compute ptrans
ptrans=pi_param*eye(nums);

ptrans(1:end-1,2:end)=ptrans(1:end-1,2:end)+((1-pi_param)/2)*eye(nums-1);
ptrans(2:end,1:end-1)=ptrans(2:end,1:end-1)+((1-pi_param)/2)*eye(nums-1);
ptrans(1,2)=1-pi_param;
ptrans(end,end-1)=1-pi_param;

%% Compute gradient
grad=eye(nums);

grad(1:end-1,2:end)=grad(1:end-1,2:end)+(-1/2)*eye(nums-1);
grad(2:end,1:end-1)=grad(2:end,1:end-1)+(-1/2)*eye(nums-1);
grad(1,2)=-1;
grad(end,end-1)=-1;

end
