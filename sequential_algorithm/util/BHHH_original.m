
%https://github.com/fediskhakov/MPECvsNFXP/blob/8212346a0426f54977823bd754947f528ed749ba/nfxp_Rust87/nfxp.m#L315

 function [mphat, mpse, cov, g, ll, iterinfo]=maxlik(llfun, p0, pnames, options);
            %NFXP.MAXLIK:  Routine to maximize likelihood functions based on
            %               - Analytical gradients.
            %               - BHHH Algorithm
            %               - Line search (step-halfing)
            %
            %   Syntax:   [mphat, mpse, cov, g, ll, convflag]=maxlik(llfun, p0, pnames, options);
            %
            %
            %   OUTPUTS:
            %       mphat:      Estimated parameters (structure with as many fields as p0)
            %       mpse:       Standard errors of estimated parameters (structure with number of fields equal to number elements in pnames, nP)
            %       cov:        Variance Covariance Matrix (nP x nP matrix)
            %       g:          Score vectors (N x nP)
            %       ll:         Likelihood contributions (N x 1 vector)
            %       convflag    Bolean (true/false) idicator of convergence
            %
            %   INPUTS:
            %       llfun:      Function handle for likelihood function
            %
            %       p0:         Structure holding starting values of parameters
            %
            %       pnames:     k dimensional cell array, that holds names of the field in the parameter structure
            %                   to be estimated. Then likelihood is only maximized with respect to the parameters in pnames
            %                   and it is thus possible to fix parameters at the values set in p0
            %
            %                   Example:
            %                       Suppose parameter structure contains both p0.delta0 and p0.delta1
            %                       If pnames={'delta0','delta1'}, then both delta0 and delta1 are estimated
            %                       If pnames={'delta0'}, then only delta0 is estimated and delta1 is fixed at starting value (set in p0)
            %
            %
            %       options:    Structure holding voluntary options for ML. If unspecified default options are used
            %                   options structure has the following fields:
            %
            %           options.title (String):
            %                   A title for the model.
            %
            %           options.output (Scalar) :
            %                   0 if no output, 1 for final output, 2 for detailed iteration output (default is 2)
            %
            %           options.maxit (Scalar):
            %                   Maximum number of iterations (default is 100)
            %
            %           options.tolgdirec (Scalar):
            %                   Tolerance for gradient'direc (default is 1e-3)
            %
            %           options.lambda0:
            %                   The basic step length (default is 1).
            % See also:
            %   nfxp.solve
            
            tt=tic;
            par.title={''};
            par.maxit=100;
            par.tolgdirec=1e-3;
            par.tol_bhhh=0.1;
            par.output=1;
            par.lambda0=1;
            par.hess=0;

            p1=p0;

            global NKIter
            iterinfo.MajorIter=0;
            iterinfo.ll=0;
            iterinfo.NKIter=0;
                        
            if nargin==4
                pfields=fieldnames(options);
                for i=1:numel(pfields);
                    par.(pfields{i})=options.(pfields{i});
                end
            end
            
            lambda=par.lambda0;   % initialize stepsize
            k = length(pnames);  % Number of parameters to estimate
            
            % -------------------------
            %  ** BHHH Optimization **
            % -------------------------
            
            if (par.hess==1) 
                strhess='Analytical Hessian';
                [ll0,ev,s0, h]=llfun(p0);
                h=-h;
            elseif (par.hess==0);  
                strhess='Outer product of scores';
                [ll0,ev,s0]=llfun(p0);
                h=s0'*s0;
            end;

            iterinfo.ll=iterinfo.ll+1;
            iterinfo.NKIter=iterinfo.NKIter+NKIter;


            g=sum(s0)';
            direc=h\g;
            tolm=g'*direc;
            if tolm<0; % convex area of likelihood, Newton-Raphson moves away from max
                h=s0'*s0;
                direc=h\g;
                tolm=g'*direc;
                if par.output>=1 
                    disp('Convex area of likelihood: Switch to BHHH')
                end
            end

            iterinfo.MajorIter=0;
            for it=0:par.maxit;
                for i=1:k;
                    np=numel(p1.(char(pnames(i))));
                    p1.(char(pnames(i)))=p0.(char(pnames(i)))+ lambda*direc(i:i+np-1);                    
                end
                p1.ev=ev;
                
                if lambda==par.lambda0;   % If previous step was accepted
                     if ((par.hess==0) | (abs(tolm) > par.tol_bhhh)); 
                        strhess='Outer product of the scores';
                        [ll,ev,s]=llfun(p1);
                        h=s'*s;
                    else
                        [ll,ev,s, h]=llfun(p1);
                        strhess='Analytical Hessian';
                        h=-h;
                    end;

                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;


                else
                    [ll, ev, s]=llfun(p1); % Do we need to recalculate s, when line searching
                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;


                end
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
                    
                    
                    % ** Plot iteration info ***
                    if par.output>1;
                        fprintf('Iteration: %d   \n', iterinfo.MajorIter);
                        fprintf('grad*direc       %10.5f \n', g'*direc);
                        fprintf('||rel grad||     %10.5f \n', max(abs(g/sum(ll0))));
                        fprintf('Log-likelihood   %10.5f \n', sum(ll));
                        fprintf('Hessian update   %10s \n', strhess);
                        fprintf('Step length      %10.5f \n', lambda);
                        fprintf('---------------------------------------------------------------------\n');
                        fprintf('    Param.       Estimates      Direction         Gradient\n');
                        i=0;
                        for iP=1:numel(pnames);
                            for j=1:numel(p0.(char(pnames(iP))));
                                i=i+1;
                                fprintf('    %-10s %10.4f    %10.4f        %10.4f\n', char(pnames(iP)), p0.(char(pnames(iP)))(j), direc(i), g(i));
                            end;
                        end;
                        fprintf('---------------------------------------------------------------------\n');

                        fprintf('\n');
                    end
                    g=sum(s)';
                    iterinfo.MajorIter=iterinfo.MajorIter+1;
                    direc=h\g;
                    tolm=g'*direc;
                    tolgrad=max(abs(g/sum(ll0)));

                    if tolm<0; % convex area of likelihood, Newton-Raphson moves away from max
                        h=s'*s;
                        direc=h\g;
                        if par.output>=0
                            disp('Convex area of likelihood: Switch to BHHH')
                        end
                    end
                    lambda=par.lambda0; % and reset stepsize
                end
                
                if tolm < par.tolgdirec;  % Stopping rule
                    break;
                end
                
            end
            
            cov=inv(h); 
            se=(sqrt(diag(cov)));
            mphat=p0;
            i=0;
            for iP=1:numel(pnames);
                for j=1:numel(p0.(char(pnames(iP))))
                    i=i+1;
                    mpse.(char(pnames(iP)))(j)=se(i);
                end
            end
            
            % ---------------------------
            %  ** Plot final output  ***
            % ---------------------------
            if it<par.maxit;
                if par.output >= 1;
                    cr=corrcov(cov);
                    fprintf('---------------------------------------------------------------------\n');
                    fprintf('***                   Convergence Achieved                        ***\n');
                    fprintf('---------------------------------------------------------------------\n');
                    disp('                         _ ');
                    disp('                         \`\ ');
                    disp('                         |= | ');
                    disp('                        /-  ;.---. ');
                    disp('                  _ __.''     (____) ');
                    disp('                   `         (_____) ');
                    disp('                   _''  ._ .'' (____) ');
                    disp('                    `        (___) ');
                    disp('                   --`''------''` ');
                    fprintf('%s \n', char(par.title));
                    fprintf('Number of iterations: %d\n', iterinfo.MajorIter);
                    fprintf('grad*direc       %10.5f \n', g'*direc);
                    fprintf('Log-likelihood   %10.5f \n', sum(ll));
                    fprintf('Step length      %10.5f \n', lambda);
                    fprintf('\n');
                    fprintf('    %-13s %13s %13s %13s\n','Param.','Estimates','s.e.','t-stat');
                    fprintf('---------------------------------------------------------------------\n');
                     i=0;
                     for iP=1:numel(pnames);
                            for j=1:numel(p0.(char(pnames(iP))))
                                i=i+1;
                                myStr = char(pnames(iP));
                                if numel(p0.(char(pnames(iP))))>1; 
                                    myStr = sprintf('%-1s(%i)',char(pnames(iP)),j);
                                end;
                                fprintf('    %-13s %13.4f %13.4f %13.4f\n', myStr, p0.(char(pnames(iP)))(j), se(i), p0.(char(pnames(iP)))(j)/se(i));
                                mpse.(char(pnames(iP)))(j)=se(i);
                            end
                    end
                    fprintf('---------------------------------------------------------------------\n');

                    fprintf('\n');

                    fprintf('Correlation matrix of parameters\n');
                    fprintf('%18s','');
                    i=0;
                    for iP=1:k;
                        for j=1:numel(p0.(char(pnames(iP))))
                            i=i+1;
                            myStr = char(pnames(iP));
                            if numel(p0.(char(pnames(iP))))>1; 
                                myStr = sprintf('%-s(%i)',char(pnames(iP)),j);
                            end;

                            fprintf('%18s', myStr);
                        end;
                    end
                    fprintf('\n');

                    i=0;
                    for iP=1:k;
                        for j=1:numel(p0.(char(pnames(iP))))
                            i=i+1;
                            myStr = char(pnames(iP));
                            if numel(p0.(char(pnames(iP))))>1; 
                                myStr = sprintf('%-1s(%i)',char(pnames(iP)),j);
                            end;

                            fprintf('%-18s', myStr);
                            for j2=1:size(cr,1);
                                fprintf('%18.4f',cr(i,j2));
                            end
                            fprintf('\n');
                        end;
                    end
                    if ismac & par.output >= 3 
                        !say 'Success! Convergence Achieved'
                    end
                    tt_end=toc(tt);
                    fprintf('\n');
                    fprintf('Time to convergence is %3.0f min and %2.2f seconds\n', floor((tt_end)/60),tt_end - 60*floor((tt_end)/60));

                    fprintf('\n');
                    
                end
                % output a flag that we did converge
                iterinfo.convflag = true;
            else
                fprintf('----------------------------------------------------------------------------\n');
                fprintf('***   BHHH failed to converge: Maximum number of iterations is reached *** \n');
                fprintf('----------------------------------------------------------------------------\n');
                disp('                   _,....._ ')
                disp('                  (___     `''-.__ ')
                disp('                 (____ ')
                disp('                 (____ ')
                disp('                 (____         ___ ')
                disp('                      `)   .-''` ')
                disp('                      /  .'' ')
                disp('                     | =| ')
                disp('                      \_\ ')
                
                if ismac
                    !say 'B triple H failed to converge: Maximum number of iterations is reached without convergence. Better luck next time'
                end
                
                % no luck this time
                iterinfo.convflag = false;
                
            end

        end % end of estim.maxlik
        
    end % end of methods

