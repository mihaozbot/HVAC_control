clear all; close all; clc;
set(0,'defaultTextInterpreter','latex')
warning('on', 'globaloptim:particleswarm:initialSwarmNotInBounds');
set(0,'defaultTextInterpreter','latex')
color = getcolors(7,"default");

HVAC_par %Load model

Q_P = Q_P_data;

c_P = 1; %People inforamtion active
c_O = 1; %Outside inforamtion active

k_0 = max([M.d_I,M.d_D])+2;

Nsim = 10000;
w = 295*ones(Nsim,1);
w(2000:5000) = 296;
w(8000:10000) = 294;

dw = 2;
wlb = w - dw;
wub = w + dw;


%  model procesa 1. reda
Model = M;
Model.N = 0; %0.5000;
Model.K_D = 0.7*M.K_D;  %11.0657;
Model.tau_D = 0.7*M.tau_D;  %1200;
% Model.K_I = M.K_I;  %1;
% Model.tau_I = M.tau_I;  %1800;
%Model.a_I = M.a_I;  %3;
Model.b_I = 0.7*M.b_I;  %-117.2899;
% Model.u_1 = M.u_1;  %-2.5000;
% Model.u_2 = M.u_2; %-1.0000e-03;
% Model.u_3 = M.u_3; %80.0000;
% Model.d_D = M.d_D; %1;
% Model.d_I = M.d_I; %4;

Q_I_m = 0*ones(k_0,1);
Q_D_m = 0*ones(k_0,1);
T_z_m = w(k_0)*ones(k_0,1);
T_O_m = c_O*T_o_data + (1-c_O)*295;
Q_P_m = c_P*Q_P;

M.act_P = 0;
Model.act_P = 0; %0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HVAC_par

k_0 = max([M.d_I,M.d_D])+2;
u = 0*ones(k_0,1);
Q_I = 0*ones(k_0,1);
Q_D = 0*ones(k_0,1);
T_z = T_z_data(1:k_0);
k_min = 1000;
sum_T_z = 0*ones(k_0,1);
Q_P = Q_P_data;
u_sum = 0*ones(k_0,1);

for k = k_0:1:(length(w))


    sum_T_z(k) = sum_T_z(k-1) + (w(k)-T_z(k-1));
    u(k) = M.u_1*(w(k)-T_z(k-1))+  M.u_2*sum_T_z(k);

    [T_z(k), Q_I(k), Q_D(k)] = HVAC(T_z(k-1), T_z(k-M.d_I-1), T_z(k-M.d_D-1), u(k-1-M.d_I), T_o_data(k-M.d_I-1), T_o_data(k-M.d_D-1), Q_I(k-1), Q_D(k-1) , Q_P(k-1),  M);
    u_sum(k) = u_sum(k-1)+u(k)^2;
end

if 1

    h1 = figure(1);
    subplot(4,1,1:2); hold off;
    p1 = plot(T_z,'b'); hold on;
    p2 = plot(wub,'k');
    plot(wlb,'k')
    p3 = plot(w,'k--');
    ax = gca;
    ax.XAxis.Exponent = 4;
    %     hcaX=ax.XRuler;
    %     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.04;
    %     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
    xlim('tight')
    ylabel('Temperature')
    legend([p1,p2,p3],'System output', 'Reference interval','Reference value','Location','best')

    subplot(4,1,3); hold off; hold off;
    h2 = plot(u,'b');
    ax = gca;
    ax.XAxis.Exponent = 4;
    %     hcaX=ax.XRuler;
    %     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.024;
    %     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;
    ylabel({'Control','signal'})
    xlim('tight')

    subplot(4,1,4); hold off;
    plot(u_sum,"color",color(1,:));
    text(length(u_sum),u_sum(end),['$\leftarrow$',num2str(ceil(u_sum(end)))],'FontSize',8)
    ax = gca;
    ax.XAxis.Exponent = 4;
    ylabel('Cost')
    legend('Sum of squared inputs','Location','best')
    xlabel('Time step')
    xlim('tight')

    pause(0.001)%
    %                 hca=gca;
    %     hcaX=hca.XRuler;
    %     hcaX.SecondaryLabel.Position(2) = hcaX.SecondaryLabel.Position(2)*1.024;
    %     hcaX.SecondaryLabel.Position(1) = hcaX.SecondaryLabel.Position(1)*1.04;

figure(1);
name = '..\HVAC_Images\HVAC_PID.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');
pause(0.0001)


end
