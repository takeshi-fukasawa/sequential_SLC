
function [param_nfxp,iter_info_nfxp]=...
    run_nfxp_func(param_initial,spec_inner_iter,spec_nfxp,varargin)


   global x_initial_cell_global count_Phi_eval

   x_initial_cell_global=[];

    n_vars=varargin{end};
    for m=1:n_vars
        x_initial0_cell{m}=varargin{m};
    end


    switch_to_BHHH_spec=1;

    count_Phi_eval=0;

    c=n_vars;
    c=c+1;opt_options=varargin{c};
    c=c+1;obj_func=varargin{c};
    c=c+1;Phi_func=varargin{c};

    c=c+1;other_input_obj_func=varargin{c};
    c=c+1;other_input_Phi_func=varargin{c};

    tStart_NFXP=tic;

    if isempty(opt_options)==1
        opt_options=optimset('GradObj','off','FinDiffType','central','MaxIter',100,...
        'Display','iter','TolFun',1e-5,...
        'TolX',1e-5,MaxFunEvals = 2000);%%%%%%%
    end

    
    if spec_nfxp.bhhh_spec==1 %Switch to BHHH later
        opt_options.TolFun=1e-1;
        opt_options.TolX=1e-1;
        %opt_options.MaxIter=1;
    end


    if isfield(spec_nfxp,'param_min')==1
        param_min=spec_nfxp.param_min;
    else
        param_min=[];
    end
    if isfield(spec_nfxp,'param_max')==1
        param_max=spec_nfxp.param_max;
    else
        param_max=[];
    end


    if isfield(spec_nfxp,'param_min')==1 | isfield(spec_nfxp,'param_max')==1

        if 1==1 % Use fmincon (not knitro)
            [param_nfxp,fval_NFXP,exitflag_NFXP,output_opt_NFXP]=...
                    fmincon(@obj_nfxp_func,param_initial,[],[],[],[],...
                            param_min,param_max,[],opt_options,...
                            Phi_func,obj_func,...
                            other_input_Phi_func,other_input_obj_func,spec_inner_iter,x_initial0_cell{:});
        else
            options = knitro_options();
            options.algorithm = 1; 
            options.outlev = 2; % ログ出力レベルを詳細に設定
            %options.hessopt=2;%BFGS
            options.gradopt=3;%central finite differenting
            %options.gradopt=1;%Exact derivative
            options.linsolver=4;
            options.bar_directinterval = 100000;
            options.blasoption = 0;

            [theta_rep_knitro,fval_knitro,exitflag_knitro,output_opt_NFXP,grad_knitro] =...
                knitro_nlp(@obj_nfxp_func,param_initial,[],[],[],[],param_min,param_max,[],[],options,[],...
                Phi_func,obj_func,...
                other_input_Phi_func,other_input_obj_func,spec_inner_iter,x_initial0_cell{:});

            param_nfxp=theta_rep_knitro;
            fval_NFXP=fval_knitro;
            exitflag_NFXP=exitflag_knitro;
        end
    else

        if 1==0 % Use fminunc (not knitro)
            tic
            [param_nfxp,fval_NFXP,exitflag_NFXP,output_opt_NFXP]=...
                    fminunc(@obj_nfxp_func,param_initial,opt_options,...
                    Phi_func,obj_func,...
                    other_input_Phi_func,other_input_obj_func,spec_inner_iter,x_initial0_cell{:});
            t_fminunc=toc;
            
        else% Use Knitro provided by Artelys

            options = knitro_options();
            options.algorithm = 1; 
            options.outlev = 2; % ログ出力レベルを詳細に設定
            %options.hessopt=2;%BFGS
            options.gradopt=3;%central finite differenting
            %options.gradopt=1;%Exact derivative
            options.linsolver=4;
            options.bar_directinterval = 100000;
            options.blasoption = 0;

            [theta_rep_knitro,fval_knitro,exitflag_knitro,output_opt_NFXP,grad_knitro] =...
                knitro_nlp(@obj_nfxp_func,param_initial,[],[],[],[],[],[],[],[],options,[],...
                Phi_func,obj_func,...
                other_input_Phi_func,other_input_obj_func,spec_inner_iter,x_initial0_cell{:});

            param_nfxp=theta_rep_knitro;
            fval_NFXP=fval_knitro;
            exitflag_NFXP=exitflag_knitro;
        end% fminunc or Knitro

    end
    iter_info_nfxp=[];
    iter_info_nfxp.exitflag_NFXP=exitflag_NFXP;
    iter_info_nfxp.output_opt_NFXP=output_opt_NFXP;
    iter_info_nfxp.fval_NFXP=fval_NFXP;

    if spec_nfxp.bhhh_spec==1 %BHHH
        %%% Use BHHH
        obj_func_bhhh=@obj_func_given_v_psteady_unobs_hetero_bhhh;%%%%%% DDC w/o unobs hetero
            
        ll_func=@(param) obj_nfxp_func_bhhh(param,Phi_func,obj_func_bhhh,...
                other_input_Phi_func,other_input_obj_func,spec_inner_iter,x_initial0_cell{:});
    
        t_start_BHHH=tic;
        param_initial=param_nfxp;%%%%%%
        [param_nfxp,mpse, cov, g, ll, output_opt_NFXP_bhhh]=...
            maxlik_BHHH_func(ll_func,param_initial,param_min,param_max);
        t_bhhh=toc(t_start_BHHH);

        iter_info_nfxp.fval_NFXP_BHHH=-sum(ll)/size(ll(:),1);

        
        iter_info_nfxp.output_opt_NFXP_bhhh=output_opt_NFXP_bhhh;

    end %BHHH

    t_nfxp=toc(tStart_NFXP);



    iter_info_nfxp.count_Phi_eval=count_Phi_eval;
    n_param=size(param_initial(:),1);
    if isempty(output_opt_NFXP)==0

        if isfield(opt_options,'GradObj')==1 && strcmp(opt_options.GradObj,'on')==1
            iter_info_nfxp.count_obj_eval=output_opt_NFXP.funcCount;
        else % Numerical derivative based 
            iter_info_nfxp.count_obj_eval=(output_opt_NFXP.funcCount)*(2*n_param+1);
        end
    end

    iter_info_nfxp.t_nfxp=t_nfxp;

end
