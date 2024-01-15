clc; clear all; close all;
set(groot,'defaultAxesFontSize',10)

HVAC_par_script_name = 'path_to_HVAC_par_script'; % Provide the correct path or name of the HVAC_par script

% swarm_size_input = 10;
% Deltau_input = 0.005;
% Maxu_input = 4;
% lambda_1_input = 100;
% lambda_2_input = 100;
% max_iter_input = 100;
% ar_input = 0.90;
% br_input = 0.10;
% H_input = 30;
% alpha_input = 0.95;

if 0

Deltau = 0.01;
lambda_1 = 100;
lambda_2 = 2;
H = 30;
ar = 0.90;
alpha = 0.95;
n_var_0 = 10;
var_0 = 1;

% Initial parameters
load  HVAC_PFC_optimal
% Extract parameters from the input vector
Deltau = par_optimized(1);
lambda_1 = par_optimized(2);
lambda_2 = par_optimized(3);
ar = par_optimized(4);
H = round(par_optimized(5));
alpha = par_optimized(6);
n_var_0 = round(par_optimized(7));
var_0 = par_optimized(8);

par_0 = [Deltau, lambda_1, lambda_2, ar,    H,   alpha, n_var_0, var_0];
lb = [0.005,     1,      0.1,        0.9    20,  0.9,   1,       1];
ub = [0.008,      1000,     10,      0.95,  30,  0.99,  10,      3];


    % options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', ...
    %     'PlotFcn', {@gaplotbestf},'PopulationSize', 20, 'UseParallel',true);
    global best_cost; best_cost = 150; % Initialize with a large value
    % [par_optimized, fval, exitflag, e] = ga(@objectiveFunction, length(par_0), [], [], [], [], lb, ub, [], options);

    options = optimoptions('particleswarm','maxIter', 100, 'Display','iter',...
        'SwarmSize',20,'UseParallel',false,'PlotFcn','pswplotbestf');%);

    options.InitialSwarmMatrix = par_0;
    [par_optimized, fval, exitflag, e] = particleswarm(@objectiveFunction, ...
        length(par_0), lb, ub, options);

    save HVAC_PFC_optimal
    
end

load  HVAC_PFC_optimization_backup
par_optimized = params;

% Extract parameters from the input vector
Deltau = par_optimized(1);
lambda_1 = par_optimized(2);
lambda_2 = par_optimized(3);
ar = par_optimized(4);
H = round(par_optimized(5));
alpha = par_optimized(6);
n_var_0 = round(par_optimized(7));
var_0 = par_optimized(8);

% Here, you need to define how the output of your function translates to a cost.
% For instance, if compute_HVAC_PFC_control returns a mean squared error, then:
[Yp, Ym, U, W, Ws, V, Z, u_sum, Q_D_p, Q_I_p, wlb, wub] = ...
    compute_HVAC_PFC_interval(...
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
text(length(u_sum),u_sum(end),['\leftarrow',num2str(ceil(u_sum(end)))],'FontSize',8)
ylabel('Cost')
legend('Sum of squared inputs','Location','best')
xlabel('Time step')
xlim('tight')
pause(0.001)

figure(1);
name = '..\HVAC_Images\HVAC_PFC_interval_control.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');
pause(0.0001)

figure(2);
name = '..\HVAC_Images\HVAC_PFC_interval_cost.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');
pause(0.0001)

function cost = objectiveFunction(params)
global best_cost;

% Snap parameters to their discrete increments
% params(1) = round(params(1) * 100) / 100;       % Deltau, increment of 0.01
% params(2) = round(params(2));
% params(3) = round(params(3));
% params(4) = round(params(4) * 100) / 100;       % ar, increment of 0.01
params(5) = round(params(5) / 1) * 1;           % H, increment of 5
% params(6) = round(params(6) * 100) / 100;       % alpha, increment of 0.01
params(7) = round(params(7));                   % n_var_0, increment of 1
% params(8) = round(params(8)*10)/ 10;            % var_0, increment of 0.1

% Extract parameters from the input vector
Deltau = params(1);
lambda_1 = params(2);
lambda_2 = params(3);
ar = params(4);
H = params(5);
alpha = params(6);
n_var_0 = params(7);
var_0 = params(8);

% Here, you need to define how the output of your function translates to a cost.
% For instance, if compute_HVAC_PFC_control returns a mean squared error, then:
[Yp, Ym, U, W, Ws, V, Z, u_sum, Q_D_p, Q_I_p, wlb, wub] = ...
    compute_HVAC_PFC_interval(...
    Deltau, lambda_1, lambda_2, ar, H, alpha, n_var_0, var_0);

% If compute_HVAC_PFC_control does not directly return a scalar value you want to minimize,
% you'll need to modify this accordingly.
cost_u2 = u_sum(end);
cost_CI = any((Yp+Z.*V)>wub(1:length(Yp))) + any((Yp-Z.*V)<wlb(1:length(Yp)));
cost = cost_u2 + 100*cost_CI;

if cost < best_cost
    best_cost = cost_u2;
    save HVAC_PFC_optimization_backup
end

end

