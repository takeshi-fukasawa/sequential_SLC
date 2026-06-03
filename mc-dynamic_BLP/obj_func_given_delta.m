function obj=obj_func_given_delta(delta,PX,X,Z,W)
    theta1 = PX*delta(:);
	resid = delta(:) - X*theta1;
	g = resid'*Z;%/size(resid(:),1);
	obj = g*W*g';
end