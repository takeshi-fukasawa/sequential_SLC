function out=construct_vars_from_phat_func(phat,...
    beta,u1_0,u2_0)

    phat1=phat(:,:,1);
    phat2=phat(:,:,2);

    nstates=size(phat1,1);

            % Note: not actually psi . . . this is e(p)
            psi1 = normpdf(norminv(phat1))./(2.*phat1);
            psi2 = normpdf(norminv(phat2))./(2.*phat2);

            %u_i = h_i*theta
            h1 = [ones(nstates,1), phat2(:,2), [1;1;0;0]];
            h2 = [ones(nstates,1), phat1(:,2), [1;0;1;0]];
            
            %conditional transition matrices
            F1_0 = [phat2(:,1), phat2(:,2), zeros(nstates,2)]; %agent 1, a = 0
            F1_1 = [zeros(nstates,2), phat2(:,1), phat2(:,2)]; %agent 1, a = 1
            F2_0 = [phat1(:,1), zeros(nstates,1), phat1(:,2), zeros(nstates,1)]; %agent 2, a = 0
            F2_1 = [zeros(nstates,1), phat1(:,1), zeros(nstates,1), phat1(:,2)]; %agent 2, a = 1
            
            %unconditional transition matrix
            F = [phat1(:,1).*phat2(:,1), phat1(:,1).*phat2(:,2), ...
               phat1(:,2).*phat2(:,1), phat1(:,2).*phat2(:,2)];
            
            %objects for NPL approximation to V
            %u[i]_1-u[i]_0=H[i]*theta+z[i]
            DG = eye(nstates)-beta.*F;
            %DGinv = inv(DG);
            H1 = h1+beta.*(F1_1-F1_0)*(DG\(phat1(:,2).*h1));
            z1 = -u1_0+beta.*(F1_1-F1_0)*(DG\(phat1(:,1).*u1_0+sum(phat1.*psi1,2)));
            H2 = h2+beta.*(F2_1-F2_0)*(DG\(phat2(:,2).*h2));
            z2 = -u2_0+beta.*(F2_1-F2_0)*(DG\(phat2(:,1).*u2_0+sum(phat2.*psi2,2)));

            out={H1,H2,z1,z2};
end
