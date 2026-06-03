%% Construct instruments, following Sun and Ishihara (2019)

x_jt_competitors=sum(reshape(x_jt,J,1,G,T,3),[1,3])-reshape(x_jt,J,1,G,T,3);%J*1*G*T*3

Z_jt_temp=[reshape(x_jt,[],3),reshape(z_jt,[],1),...
    reshape(w_jt,[],1)];%(J*G*T)*9
Z_jt_temp=[Z_jt_temp,reshape(x_jt_competitors,[],3)];

Z_jt_temp=Z_jt_temp./max(abs(Z_jt_temp));


Z_jt_temp_quadratic=NaN(J*G*T,1);
Z_jt_temp_cubic=NaN(J*G*T,1);

id=1;
for m=1:size(Z_jt_temp,2)
   for n=1:m
       Z_jt_temp_quadratic(:,id)=Z_jt_temp(:,m).*Z_jt_temp(:,n);
       id=id+1;
    end % loop wrt n
end % loop wrt m

id=1;
for m=1:size(Z_jt_temp,2)
   for n=1:m
       for q=1:n
           Z_jt_temp_cubic(:,id)=Z_jt_temp(:,m).*Z_jt_temp(:,n).*Z_jt_temp(:,q);
           id=id+1;
       end
    end % loop wrt n
end % loop wrt m

%%% Differentiation IV
K_x=size(x_jt,2);
d_mat=reshape(x_jt,J*G,1,T,K_x)-reshape(x_jt,1,J*G,T,K_x);%(J*G)*(J*G)*T*K_x
temp=sum(d_mat.^2,2);%(J*G)*1*T*K_x
diff_IV1=reshape(temp,J*G*T,K_x);

diff_IV=[];
for k1=1:K_x
    for k2=1:k1
        temp=sum(d_mat(:,:,:,k1).*d_mat(:,:,:,k2),2);%(J*G)*1*T
        diff_IV=cat(2,diff_IV,reshape(temp,J*G*T,1));
    end
end

diff_IV=diff_IV1;%%%%
diff_IV=diff_IV./max(abs(diff_IV));

w_jt=reshape(w_jt,[],1);
Z_jt=[ones(J*G*T,1),Z_jt_temp];%(J*G*T)*n_IV
%Z_jt=[Z_jt,Z_jt_temp_quadratic,Z_jt_temp_cubic];
%Z_jt=[Z_jt,diff_IV,w_jt.^2,w_jt.*x_jt];
%Z_jt=[Z_jt,diff_IV];%,Z_jt_temp_quadratic];
Z_jt=[Z_jt,diff_IV,Z_jt_temp_quadratic];
%Z_jt=[Z_jt,Z_jt_temp_quadratic,Z_jt_temp_cubic];

%Z_jt=Z_jt(:,1:5+3);

X=reshape(X_jt,[],size(X_jt,5));
Z=Z_jt;

W=inv(Z'*Z);

PX = inv(X'*Z*W*Z'*X)*X'*Z*W*Z';
