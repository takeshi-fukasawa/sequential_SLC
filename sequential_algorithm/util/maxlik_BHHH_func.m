
%https://github.com/fediskhakov/MPECvsNFXP/blob/8212346a0426f54977823bd754947f528ed749ba/nfxp_Rust87/nfxp.m#L315: 
%%% Replication code for Iskhakov et al. (2015) "Constrained Optimization Approaches To Estimation Of Structural Models: Comment" in Econometrica
%%% Modified by Takeshi Fukasawa in Feb 2026 (Introduce parameter bounds)

 function [mphat, mpse, cov, g, ll, iterinfo]=maxlik_BHHH_func(llfun, p0,param_min,param_max);
            %maxlik:  Routine to maximize likelihood functions based on
            %               - Analytical gradients.
            %               - BHHH Algorithm
            %               - Line search (step-halfing)
            %
            %
            %   OUTPUTS:
            %       mphat:      Estimated parameters
            %       mpse:       Standard errors of estimated parameters             %       cov:        Variance Covariance Matrix (nP x nP matrix)
            %       g:          Score vectors (N x nP)
            %       ll:         Likelihood contributions (N x 1 vector)
            %       convflag    Bolean (true/false) idicator of convergence
            %
            %   INPUTS:
            %       llfun:      Function handle for likelihood function
            %
            %       p0:        Starting values of parameters
            %
            %            
            tt=tic;
            par.title={''};
            par.maxit=100;
            par.tolgdirec=1e-5;%1e-3;
            par.tol_bhhh=0.1;
            par.output=2;
            par.lambda0=1;

            p1=p0;

            global NKIter
            iterinfo.MajorIter=0;
            iterinfo.ll=0;
            iterinfo.NKIter=0;
                        
            lambda=par.lambda0;   % initialize stepsize
            
            % -------------------------
            %  ** BHHH Optimization **
            % -------------------------
            
            [ll0,s0]=llfun(p0);
            h=s0'*s0;
            n_obs=size(ll0(:),1);

            iterinfo.ll=iterinfo.ll+1;
            iterinfo.NKIter=iterinfo.NKIter+NKIter;


            g=sum(s0)';
            direc=h\g;
            %if norm(direc)>=2
            %    direc=direc/norm(direc);
            %end

            tolm=g'*direc;

            iterinfo.MajorIter=0;
            for it=0:par.maxit;
                p1=p0+lambda*direc;
                if isempty(param_min)==0 
                    p1=max(p1,param_min);
                end
                if isempty(param_max)==0 
                    p1=min(p1,param_max);
                end
                
                [ll,s]=llfun(p1);
                h=s'*s;
                iterinfo.ll=iterinfo.ll+1;
                iterinfo.NKIter=iterinfo.NKIter+NKIter;

                if sum(ll)<sum(ll0);
                    lambda=lambda/2;
                    if par.output>1;
                        fprintf('\n');
                        fprintf('.............LINESEARCHING\n');
                        fprintf('\n');
                    end
                    if lambda<=.01;
                        if par.output>1;
                            fprintf('.............WARNING: Failed to increase - downhill step accepted\n');
                        end
                        p0=p1;
                        ll0=ll;                 % Accept step
                        lambda=par.lambda0;
                    end
                else       % IF INCREASE IN LIKELIHOOD, CALCULATE NEW DIRECTION VECTOR
                    p0=p1;
                    ll0=ll;                 % Accept step
                    
                    tolgrad=max(abs(g/(sum(ll0)/n_obs)));
                    
                    % ** Plot iteration info ***
                    if par.output>1;
                        fprintf('Iteration: %d   \n', iterinfo.MajorIter);
                        fprintf('grad*direc       %10.5f \n', g'*direc);
                        fprintf('||rel grad||     %10.5f \n', tolgrad);
                        fprintf('Log-likelihood   %10.5f \n', sum(ll));
                        fprintf('Step length      %10.5f \n', lambda);
                        fprintf('\n');
                    end
                    g=sum(s)';
                    iterinfo.MajorIter=iterinfo.MajorIter+1;
                    direc=h\g;
                    tolm=g'*direc;

                    %tolgrad=max(abs(g));
                    lambda=par.lambda0; % and reset stepsize
                end
                
                %if tolm < par.tolgdirec;  % Stopping rule
                if tolgrad < par.tolgdirec | max(abs(direc))<par.tolgdirec;  % Stopping rule
                    break;
                end
                
            end%for it=0:par.maxit;
            
            cov=inv(h); 
            se=(sqrt(diag(cov)));
            mphat=p0;
            mpse=se;
    end 
        
