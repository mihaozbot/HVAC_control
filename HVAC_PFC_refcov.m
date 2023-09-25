clc; clear all; close all;

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

%  proces 1. red
Q_I_p = 0*ones(k_0,1);
Q_D_p = 0*ones(k_0,1);
T_z_p = w(k_0)*ones(k_0,1);
T_O_p = T_o_data;
Q_P_p  = Q_P;

%  model procesa 1. reda
Model = M;
Model.N = 0; %0.5000;
Model.K_D = 0.8*M.K_D;  %11.0657;
Model.tau_D = 0.8*M.tau_D;  %1200;
% Model.K_I = M.K_I;  %1;
% Model.tau_I = M.tau_I;  %1800;
%Model.a_I = M.a_I;  %3;
Model.b_I = 0.8*M.b_I;  %-117.2899;
% Model.u_1 = M.u_1;  %-2.5000;
% Model.u_2 = M.u_2; %-1.0000e-03;
% Model.u_3 = M.u_3; %80.0000;
% Model.d_D = M.d_D; %1;
% Model.d_I = M.d_I; %4;

Q_I_m = 0*ones(k_0,1);
Q_D_m = 0*ones(k_0,1);
T_z_m = w(k_0)*ones(k_0,1);
T_O_m = c_O*T_o_data + (1-c_O)*294;
Q_P_m = c_P*Q_P;

M.act_P = 1;
Model.act_P = 1; %0;

% parametri simulacije
yp = T_z_p(end); ym = T_z_m(end);

z_sig = 5;
n_var_0 = 10;
n_var = n_var_0;
sigma_0 = n_var*((dw/z_sig)/2)^2;
var2 = sigma_0;
sig = sqrt(var2/n_var);
mu = 0; %(wlb(1) + 3*sqrt(sigma2/n_var));

U = 0*ones(k_0,1); Ym = T_z_m; Yp= T_z_p;  W = T_z_p;  uZ=0; Ws = T_z_p;
V = sig*ones(k_0,1);

% referenèni model
ar = 0.90; br = (1-ar); H = 10;
N1 = 1;
N2 = H;
u_sum = 0*ones(k_0,1);
% ---- inicializacija -----
swarm_size = 20;
for i = 1 : swarm_size
    swarm(i, 1, 1) = i;
end

u = 0; Maxu = 6; Deltau = 0.5;

t1 = clock;

for iSim = k_0:Nsim

    w_sigma = w(iSim); %min((wlb(iSim) + z_sig*sig),w(iSim));
    %e = w(iSim)-yp;

    %ymRefH2=ym+e*(1-ar^H/2);
    % PSO algoritem
    iterations = 100;
    inertia = 1.0;
    correction_factor = 2;

    % ---- inicializacija -----
    for i = 1 : swarm_size
        %      swarm(i, 1, 1) =  u - round(swarm_size/2)*Maxu/swarm_size + i * (Deltau);
        swarm(i, 1, 1) =  i * (Maxu/swarm_size)*0.1 + u ;
    end

    %     for i = 1 : swarm_size
    %         %      swarm(i, 1, 1) =  u - round(swarm_size/2)*Maxu/swarm_size + i * (Deltau);
    %         swarm(i, 1, 2) =  i * ((wub(iSim)- wlb(iSim))/swarm_size)*0.1 + w(iSim);
    %     end

    %     swarm_dw = repelem(( min(wlb(iSim) + z_sig*sig, w(iSim))-w(iSim)):...
    %         ((max(wub(iSim) - z_sig*sig,w(iSim)) - min(wlb(iSim) + z_sig*sig, w(iSim)))/(swarm_size/2-1)):...
    %         (max(wub(iSim) - z_sig*sig,w(iSim))-w(iSim)),1,2);
    %
    % %     if ~isempty(swarm_dw)
    % %
    % %
    % %     else
    % %         swarm(:, 1, 2) = w(iSim);
    % %     end
    %

    %      swarm(i, 1, 1) =  u - round(swarm_size/2)*Maxu/swarm_size + i * (Deltau);
    swarm(:, 1, 2) =  repmat(floor( 1:swarm_size-swarm_size/2).*((wub(iSim)-wlb(iSim))/swarm_size)*1 + w(iSim),1,2);
    swarm((swarm(:, 1, 2) < min(wlb(iSim) + z_sig*sig, w(iSim))) |...
        (swarm(:, 1, 2)  > max(wub(iSim) - z_sig*sig,w(iSim))), 1, 2) = w(iSim);

    swarm(:,3,1) = swarm(:,1,1);
    swarm(:,3,2) = swarm(:,1,2);

    swarm(:, 4, 1) = 10000;          % najboljša vrednsot do sedaj
    swarm(:, 4, 2) = NaN;          % najboljša vrednsot do sedaj
    swarm(:, 2, :) = 0;             % zaèetna hitrost

    term = 1000; termMax = 0.01;
    % izdelaj termination criteria
    lambda = 0.1;
    gamma = 100;
    max_iter = 1000;
    iter = 0;
    while ((term > termMax) && (iter<max_iter))
        %for iter = 1 : iterations
        iter = iter +1;
        %-- evaluacija pozicije in kriterijske funkcije ---
        for i = 1 : swarm_size

            swarm(i, 1, :) = swarm(i, 1, :) + swarm(i, 2,:);     % popravek pozicije x
            ux = swarm(i, 1, 1);
            wx = swarm(i, 1, 2);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h0 = max([M.d_I,M.d_D])+2-1;
            T_z_m_H(1:h0) = T_z_m((end-h0+1):end);
            Q_I_m_H = Q_I_m((end-h0+1):end);
            Q_D_m_H = Q_D_m((end-h0+1):end);
            T_O_m_H = T_O_m((iSim-h0+1):iSim+H-1);
            Q_P_m_H = Q_P_m((iSim-h0+1):iSim+H-1);
            for h = (h0+1):1:(h0+H)
                [T_z_m_H(h), Q_I_m_H(h), Q_D_m_H(h)] = HVAC(T_z_m_H(h-1), T_z_m_H(h-Model.d_I-1),  T_z_m_H(h-Model.d_D-1), ux, ...
                    T_O_m_H(h-Model.d_I-1), T_O_m_H(h-Model.d_D-1), Q_I_m_H(h-1), Q_D_m_H(h-1), Q_P_m_H(h-1), Model);
            end
            ymH = T_z_m_H(h0+H);
            %ymH2 = T_z_m_H(h0+H/2);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ePSO = wx - yp;
            em = (yp - T_z_m_H(h0)) + T_z_m_H;
            ymRefH=ym+ePSO*(1-ar^H);

            if abs(ux) > Maxu
                weightU = 10000;
            else
                weightU = 0;
            end
            %
            if ((wx < min(wlb(iSim) + z_sig*sig, w(iSim))) || (wx > max(wub(iSim) - z_sig*sig,w(iSim))))  %abs(w-wx) > dw
                weightW = 10000;
            else
                weightW = 0;
            end

            %
            elb = wlb(iSim) -(em - z_sig*sig);
            elb(elb<0) = 0;
            eub = (em + z_sig*sig) - wub(iSim);
            eub(eub<0) = 0;
            eb = elb + eub;

            %             if any(((em - z_sig*sig) < wlb(iSim)) | ((em + z_sig*sig) > wub(iSim)))
            %                 weightYm = 100;
            %             else
            %                 weightYm = 0;
            %             end

            if abs(ux-uZ) > Deltau
                weightUz = 10000;
            else
                weightUz = 0;
            end

            val = (ymH-ymRefH).^2  + lambda*sum(ux.^2) + gamma*sum(eb.^2) + weightU + weightUz + weightW; % izraèun kriterijske funkcije
            %val= sum(((T_z_m_H(h0+N1:h0+N2-1)-T_z_m_H(h0+N1+1:h0+N2))-(ym+e*(1-ar.^(N1:N2-1))-(ym+e*(1-ar.^(N1+1:N2))))).^2) + weightU + weightUz;

            if val < swarm(i, 4, 1)                 % èe je nova pozicija boljša
                swarm(i, 3, 1) = swarm(i, 1, 1);    % popravek pozicije x,
                swarm(i, 3, 2) = swarm(i, 1, 2);    % popravek pozicije x,
                swarm(i, 4, 1) = val;
                T_z_m_H_best = T_z_m_H; % najboljša vrednsot delca - pBest
            end
        end

        swarm(1:swarm_size,1,1);
        term = std(swarm(:,1,2)) + std(swarm(:,1,1));

        [temp, gbest] = min(swarm(:, 4, 1));        % gBest

        %-- popravek pozicije delca
        for i = 1 : swarm_size
            swarm(i, 2, :) = rand*inertia*swarm(i, 2, :) + correction_factor*rand*(swarm(i, 3, :) - swarm(i, 1, :)) ...
                + correction_factor*rand*(swarm(gbest, 3, :) - swarm(i, 1, :));   %x komponenta hitrosti
        end
        %
        %         %         %% Plotanje swarma
        %                     clf
        %                     plot(swarm(:, 1, 1),swarm(:, 1, 2)-w(iSim), 'x')   %  izris premikanja delca%
        % %                     xlim([-5,5]);
        % %                     ylim([-5,5]);
        % %                      pause(.2)
    end

    u = swarm(gbest,1,1);
    uZ = u;

    wZ = swarm(gbest,1,2);

    e = wZ - yp;
% 
%     if (wZ + z_sig*sig) < yp
%         sig = (yp - wZ)/(z_sig)/2;
%         var2 = (sig^2)*n_var;
%     end

    if 0
        figure(2); hold off;
        plot(0:h0,Yp((end-h0):end),'b'); hold on;
        plot(h0:h0+H, em(h0:h0+H),'b--')
        plot(0:h0,Ym((end-h0):end),'r')
        plot(h0:h0+H,T_z_m_H_best(h0:h0+H),'r--');
        plot(h0:h0+H,ym+e*(1-ar.^(0:H)),'k');
        plot(h0:h0+H,[ym, repelem(ymRefH,1,H)],'m--');
        xline(h0)
        legend('Process past','Estimated process future','Model past','Model future','Reference future','Reference PFC')
        pause(0.0001)
    end

    % simulacija procesa
    d = -0.2; noise = M.N;
    up = u;
    u_sum(iSim) = u_sum(iSim-1)+up^2;
    %     if iSim > round(2*Nsim/4)
    %         up = u + d;
    %     end

    %yp = ap*ypZ + bp*up;
    %yp = yp + noise*randn(1,1);
    %ypZ = yp;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [T_z_p(iSim), Q_I_p(iSim), Q_D_p(iSim)] = HVAC(T_z_p(iSim-1), T_z_p(iSim-M.d_I-1), T_z_p(iSim-M.d_D-1), up,...
        T_O_p(iSim-M.d_I-1), T_O_p(iSim-M.d_D-1), Q_I_p(iSim-1), Q_D_p(iSim-1), Q_P_p(iSim-1), M);
    yp = T_z_p(iSim) + noise*randn(1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if abs(w(iSim)- w(iSim-1))>0
        %w_sigma = w(iSim);
        var2 = sigma_0;
        n_var = n_var_0;
        mu = e;
        T_z_m((end-h0):end) = T_z_p((end-h0):end);

    end

    e_var = e - mu; %Distance of new data from center
    mu = mu + 1/(1 + n_var)*e_var; %Center update
    var2  = var2 + e_var*(e - mu)'; %un-normalized covariance matrix
    n_var = n_var + 1; %Increase number of samples in cluster
    sig = sqrt(var2/n_var);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % simulacije modela
    um = u;
    % ym=am*ymZ+bm*um;
    % ymZ=ym;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [T_z_m(iSim), Q_I_m(iSim), Q_D_m(iSim)] = HVAC(T_z_m(iSim-1), T_z_m(iSim-Model.d_I-1), T_z_m(iSim-Model.d_D-1), um, ...
        T_O_m(iSim-Model.d_I-1), T_O_m(iSim-Model.d_D-1), Q_I_m(iSim-1), Q_D_m(iSim-1), Q_P_m(iSim-1), Model);
    ym = T_z_m(iSim);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Yp = [Yp; yp];
    Ym = [Ym; ym];
    U = [U; up];
    W = [W; w(iSim)];
    Ws = [Ws; wZ];
    V = [V; sig];
    if (mod(iSim,100) == 0)

        h1 = figure(1);
        set(h1,'Position',[680,558,2*560,420])
        subplot(3,1,1:2); hold off;
        p1 = plot(Yp,'b'); hold on;
        p2 = plot(Ws + z_sig*V,'-.',"color",color(1,:));
        plot(Ws - z_sig*V,'-.',"color",color(1,:));
        p3 = plot(Ym,'r');
        p4 = plot(wub,'k');
        plot(wlb,'k')
        %plot(w + 2*dw,'r--')
        %plot(w - 2*dw,'r--')
        %plot(Ws,'k');


        legend([p1,p2,p3,p4],'System output','Confidence interval','Model output','Reference interval','')
        xlim('tight')
        xlabel('Time step')


        h2 = figure(2);
        subplot(2,1,1); hold off;
        h2 = plot(U,'b'); hold on;

        plot(normalize(Q_D_p(1:length(U)-1) + M.act_P*Q_P_p(1:length(U)-1)'),'r');
        legend('Control signal', 'Disturbances energy')
        xlim('tight')
        subplot(2,1,2); hold off;
        plot(u_sum,'r');
        legend('Sum of squared inputs')
        xlabel('Time step')
        pause(0.001)
        xlim('tight')

    end
end
t2 = clock;

Povprecni_Cas_Enega_Izracuna = etime(t2,t1)/Nsim

figure(1);
name = '..\HVAC_Images\HVAC_PFC_refcov_control.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');
pause(0.0001)

figure(2);
name = '..\HVAC_Images\HVAC_PFC_refcov_cost.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');
pause(0.0001)