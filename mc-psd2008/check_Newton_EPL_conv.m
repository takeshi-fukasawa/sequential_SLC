% Estimates the Pesendorfer and Schmidt-Dengler (2008) model with
% several sample sizes and values of c.
clear all
rng(0,'twister');

nsims = 100;

%for eqm_dgp = [ 1, 2, 3 ]
for eqm_dgp=2:2
    %diary(sprintf('psd2008_estimate_eqm%d.log', eqm_dgp));
    %for nobs = [ 250, 1000 ]
    for nobs=250:250
        % Generate data from the model
        generate_data

        %% EPL iteration
        EPL_not_update_param_spec=0;
        initial_est_spec=1;
        % Estimate the model with consistent CCPs
        c = 1.0;
        nstart = 1;
        psd2008_estimate
        theta_EPL_init=theta_EPL;
        save("theta_EPL_temp.mat",'theta_EPL_init')
    
        initial_est_spec=0;
        if 1==1

            %% EPL iteration
            EPL_not_update_param_spec=0;
            % Estimate the model with consistent CCPs
            c = 1.0;
            nstart = 1;
            psd2008_estimate
    
            G_val_results_EPL=G_val_results;
            save("EPL_results.mat",'G_val_results_EPL','theta_EPL_0step','theta_EPL_1step',...
                'theta_EPL_2step','theta_EPL_3step','diverge_EPL')
            
        end        

        clear theta_EPL
        if 1==1
    
            %% Newton iteration
            EPL_not_update_param_spec=1;
            % Estimate the model with consistent CCPs
            c = 1.0;
            nstart = 1;
            psd2008_estimate
    
            G_val_results_Newton=G_val_results;
            diverge_Newton=diverge_EPL;
            save("Newton_results.mat",'G_val_results_Newton','diverge_Newton')

        end

    end

    %diary off
end

load("EPL_results.mat")
diverge_newton_mat=cumsum(isnan(G_val_results_Newton));
diverge_EPL_mat=cumsum(isnan(G_val_results_EPL));

diverge_id=find(isnan(sum(G_val_results_Newton)))
nondiverge_id=find(isnan(sum(G_val_results_Newton))==0)

G_val_results_Newton_diverged=G_val_results_Newton(:,diverge_id);
G_val_results_EPL_diverged=G_val_results_EPL(:,diverge_id);
G_val_results_Newton_diverged(isnan(G_val_results_Newton_diverged))=NaN;
mean_G_newton=mean(G_val_results_Newton_diverged,2,'omitnan');
mean_G_EPL=mean(G_val_results_EPL_diverged,2,'omitnan');
std_G_newton=std(G_val_results_Newton_diverged,0,2,'omitnan');
std_G_EPL=std(G_val_results_EPL_diverged,0,2,'omitnan');

G_val_results_Newton_diverged(G_val_results_Newton_diverged==0)=10000;
plot(max(G_val_results_EPL_diverged,[],2))
plot(max(G_val_results_Newton_diverged,[],2))

temp=[mean_G_newton(1:4)',mean_G_EPL(2:4)';...
    std_G_newton(1:4)',std_G_EPL(2:4)']

diverge_rate_EPL=mean(diverge_EPL_mat,2)*100;
diverge_rate_newton=mean(diverge_newton_mat,2)*100;

temp=[diverge_rate_newton(1:9)',diverge_rate_newton(end);...
    mean_G_newton(1:9)',NaN;...
    std_G_newton(1:9)',NaN;...
 diverge_rate_EPL(1:9)',diverge_rate_EPL(end);...
    mean_G_EPL(1:9)',NaN;...
    std_G_EPL(1:9)',NaN;...
    ];
Newton_EPL_convergence_comparison_table=round(temp,2);
writematrix(Newton_EPL_convergence_comparison_table,"Newton_EPL_convergence_comparison_table.csv")


%%%%%%%%%%%%%%%
load("EPL_results.mat")
load("Newton_results.mat")

diverge_newton_mat=cumsum(isnan(G_val_results_Newton));
diverge_EPL_mat=cumsum(isnan(G_val_results_EPL));


diverge_rate_EPL=mean(diverge_EPL_mat,2)*100;
diverge_rate_newton=mean(diverge_newton_mat,2)*100;

temp=sum(G_val_results_Newton);
diverge_id=find(isnan(temp)==1);
nondiverge_id=find(isnan(temp)==0);

for i=1:size(diverge_id,2)
    id=diverge_id(i);
    G_val_results_Newton(G_val_results_Newton(:,id)==0,id)=1000;
    G_val_results_EPL(G_val_results_EPL(:,id)==0,id)=1e-6;
end
for i=1:size(nondiverge_id,2)
    id=nondiverge_id(i);
    G_val_results_Newton(G_val_results_Newton(:,id)==0,id)=1e-6;
    G_val_results_EPL(G_val_results_EPL(:,id)==0,id)=1e-6;    
end

%% Plot the convergence of EPL and Newton
fontsize=16

%set(gcf,'FontSize',fontsize)

clear title
t = tiledlayout(1,2);
nexttile
yyaxis left
temp=log10(G_val_results_EPL');
boxchart(temp(:,1:20),'MarkerStyle','.')
xlabel('iteration','FontSize',fontsize)
ylim([-6 3])
ylabel('$\log_{10}\left(\left\Vert G\left(Y;\theta\right)\right\Vert _{\infty}\right)$',...
    'Interpreter','latex','FontSize',fontsize)

    
yyaxis right
plot(diverge_rate_EPL(1:20))
ylim([0,100])
ylabel('Divergence rate (%)','FontSize',fontsize)
title('(a). SLC/EPL','FontSize',fontsize)


nexttile
yyaxis left
temp=log10(G_val_results_Newton');
p=boxchart(temp(:,1:20),'MarkerStyle','.')
xlabel('iteration','FontSize',fontsize)
ylim([-6 3])
ylabel('$\log_{10}\left(\left\Vert G\left(Y;\theta\right)\right\Vert _{\infty}\right)$',...
    'Interpreter','latex','FontSize',fontsize)

yyaxis right
plot(diverge_rate_newton(1:20))
ylim([0,100])
ylabel('Divergence rate (%)','FontSize',fontsize)
title('(b). Newton','FontSize',fontsize)

%set(gcf,'FontSize',fontsize)


%fontsize(16,"points")

saveas(gcf,'Newton_EPL_convergence.png')
