clc; clear all; close all;
set(groot,'defaultAxesFontSize',10)

HVAC_par_script_name = 'path_to_HVAC_par_script'; % Provide the correct path or name of the HVAC_par script

if 1 

global ar n_var_0 var_0 Deltau lambda_2;
Deltau = 0.005;
lambda_1 = 2;
lambda_2 = 1;
H = 100;
ar = 0.90;
alpha = 0.9;
n_var_0 = 1;
var_0 = 0;

par_0 = [lambda_1, alpha,  H];
lb = [1e-3,            0.7    30];
ub = [1e3,        0.99 , 300];

% options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', ...
%     'PlotFcn', {@gaplotbestf},'PopulationSize', 20, 'UseParallel',true);
global best_cost; best_cost = 1000; % Initialize with a large value
% [par_optimized, fval, exitflag, e] = ga(@objectiveFunction, length(par_0), [], [], [], [], lb, ub, [], options);

options = optimoptions('particleswarm','maxIter', 20, 'Display','iter',...
    'SwarmSize',100,'UseParallel',false,'PlotFcn','pswplotbestf');%);

options.InitialSwarmMatrix = par_0;
[par_optimized, fval, exitflag, e] = particleswarm(@objectiveFunction, ...
    length(par_0), lb, ub, options);

save HVAC_PFC_optimal_outside
end

load  HVAC_PFC_optimization_outside_backup
par_optimized = params;

% Extract parameters from the input vector
lambda_1 = par_optimized(1);
alpha = par_optimized(2);

% Here, you need to define how the output of your function translates to a cost.
% For instance, if compute_HVAC_PFC_control returns a mean squared error, then:
[Yp, Ym, U, W, V, Z, u_sum, Q_D_p, Q_I_p, wlb, wub] = ...
    compute_HVAC_PFC_outside(...
    Deltau, lambda_1, lambda_2, ar, H, alpha, n_var_0, var_0);

h0 = 3;
k_tics = 0:floor((h0+H)/10):(h0+H);
labels = arrayfun(@(x) sprintf('k%s', num2str(x-h0, '%+d')), k_tics, 'UniformOutput', false);
labels{find(k_tics== h0)} = 'k';
color = getcolors(7,"default");
w_shift = 200;
wlb_cor = min([[wlb(1:w_shift);wlb(1:end-w_shift)],...
    [wlb(w_shift+1:end);wlb((end-w_shift+1):end)]],[],2);
wub_cor = max([[wub(1:w_shift);wub(1:end-w_shift)],...
    [wub(w_shift+1:end);wub((end-w_shift+1):end)]],[],2);

h1 = figure(1);
set(h1,'Position',[680, 558, 2*560, 420])
subplot(3,1,1:2); hold off;
p1 = plot(Yp,'b'); hold on;
p2 = plot(Yp + Z.*V,'-.',"color",color(1,:));
plot(Yp - Z.*V,'-.',"color",color(1,:));
p3 = plot(Ym,'r');
p4 = plot(wlb_cor,'k--');
plot(wub_cor,'k--')
%p5 = plot(w,'k');
ylim('tight')
ylabel('Temperature')
legend([p1,p2,p3,p4],'System output','Confidence interval','Model output','Reference interval')
xlim('tight')
xlabel('Time step')

h2 = figure(2);
subplot(3,1,1); hold off;
h2 = plot(U,'b');
ylabel({'Control','signal'})
xlim('tight')
subplot(3,1,2); hold off;
plot((Q_D_p(1:length(U)-1)),"color",color(4,:)); hold on;
plot((Q_I_p(1:length(U)-1)),'r');
legend('Disturbance', 'Control','location','best')
ylabel({'Thermal','energy'})
xlim('tight')
subplot(3,1,3); hold off;
plot(u_sum,"color",color(1,:));hold on;
text(length(u_sum),u_sum(end),['$\leftarrow$',num2str(ceil(u_sum(end)))],'FontSize',8)
ylabel('Cost')
legend('Sum of squared inputs','Location','best')
xlabel('Time step')
xlim('tight')
pause(0.001)

function cost = objectiveFunction(params)
global ar n_var_0 var_0 Deltau lambda_2
global best_cost;

% Extract parameters from the input vector
lambda_1 = params(1);
alpha = params(2);
H = round(params(3));

% Here, you need to define how the output of your function translates to a cost.
% For instance, if compute_HVAC_PFC_control returns a mean squared error, then:
[Yp, Ym, U, W, V, Z, u_sum, Q_D_p, Q_I_p, wlb, wub] = ...
    compute_HVAC_PFC_outside(...
    Deltau, lambda_1, lambda_2, ar, H, alpha, n_var_0, var_0);

% If compute_HVAC_PFC_control does not directly return a scalar value you want to minimize,
% you'll need to modify this accordingly.
cost_u2 = u_sum(end);
cost_CI = any((Yp)>wub(1:length(Yp))) + any((Yp)<wlb(1:length(Yp)));
cost = cost_u2 + 100*cost_CI;

if cost < best_cost
    best_cost = cost_u2;
    save HVAC_PFC_optimization_outside_backup
end

%HVAC_PFC_save_all

end

