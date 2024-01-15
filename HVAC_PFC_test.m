clc; clear all; close all;
set(groot,'defaultAxesFontSize',10)
color = getcolors(7, "default");

if 0

    HVAC_par %Load model
    Q_P = Q_P_data;

    c_P = 0; %People inforamtion active
    c_O = 0; %Outside inforamtion active

    k_0 = max([M.d_I,M.d_D])+2;

    Nsim = 2000;
    w = 295*ones(Nsim,1);
    w(400:1000) = 296;
    w(1600:2000) = 294;
    dw = 2;
    wlb = w - dw;
    wub = w + dw;

    %  proces 1. red
    Q_I_p = 0*ones(k_0,1);
    Q_D_p = 0*ones(k_0,1);
    T_z_p = (w(k_0)+1)*ones(k_0,1);
    T_O_p = c_O*T_o_data + (1-c_O)*295;
    Q_P_p  = c_P*Q_P;

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

    % parametri simulacije
    yp = T_z_p(end); ym = T_z_m(end);

    z_sig = 2;
    n_var_0 = 1;
    n_var = n_var_0;
    sigma_0 = n_var*((dw/z_sig))^2;
    var2 = sigma_0;
    sig = sqrt(var2/n_var);
    mu = yp; %(wlb(1) + 3*sqrt(sigma2/n_var));

    U = 0*ones(k_0,1); Ym = T_z_m; Yp= T_z_p;  W = T_z_p;  uZ=0; Ws = T_z_p;
    V = sig*ones(k_0,1);
    Z = z_sig*ones(k_0,1);

    % referenèni model
    ar = 0.90; br = (1-ar); H = 30;
    N1 = 1;
    N2 = H;
    u_sum = 0*ones(k_0,1);

    % ---- inicializacija -----
    swarm_size = 20;
    for i = 1 : swarm_size
        swarm(i, 1, 1) = i;
    end

    u = 0; Maxu = 1; Deltau = 0.001;
    
    %name_test = 'HVAC_PFC_test';
    name_test = ['HVAC_PFC_test_delta_U_',num2str(Deltau*1e3)];

    termMax = 0.001;
    lambda_1 = 0;
    lambda_2 = 1;
    max_iter = 100;
    t1 = clock;

    h0 = max([M.d_I,M.d_D])+2-1;
    k_tics = 0:floor((h0+H)/10):(h0+H);
    labels = arrayfun(@(x) sprintf('k%s', num2str(x-h0, '%+d')), k_tics, 'UniformOutput', false);
    labels{find(k_tics== h0)} = 'k';

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
            swarm(i, 1, 1) = Deltau/2 - i * (Deltau/swarm_size) + u ;
        end

        swarm(:, 1, 2) =  repmat(floor( 1:swarm_size-swarm_size/2).*((wub(iSim)-wlb(iSim))/swarm_size)*1 + w(iSim),1,2);
        swarm((swarm(:, 1, 2) < min(wlb(iSim) + z_sig*sig, w(iSim))) |...
            (swarm(:, 1, 2)  > max(wub(iSim) - z_sig*sig,w(iSim))), 1, 2) = w(iSim);

        swarm(:,3,1) = swarm(:,1,1);
        swarm(:,3,2) = swarm(:,1,2);
        val_0 = 1e5;
        swarm(:, 4, 1) = 1e5;          % najboljša vrednsot do sedaj
        swarm(:, 4, 2) = NaN;          % najboljša vrednsot do sedaj
        swarm(:, 2, :) = 0;             % zaèetna hitrost

        % izdelaj termination criteria

        iter = 0;
        term = 1000;
        while ((term > termMax) && (iter<max_iter))
            %for iter = 1 : iterations
            iter = iter +1;
            %-- evaluacija pozicije in kriterijske funkcije ---
            for i = 1 : swarm_size

                swarm(i, 1, :) = swarm(i, 1, :) + swarm(i, 2,:);     % popravek pozicije x
                ux = swarm(i, 1, 1);
                wx = w(iSim);%swarm(i, 1, 2);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                h0 = max([M.d_I,M.d_D])+2-1;
                T_z_m_H(1:h0) = T_z_m((end-h0+1):end);
                Q_I_m_H = Q_I_m((end-h0+1):end);
                Q_D_m_H = Q_D_m((end-h0+1):end);
                T_O_m_H = T_O_m((iSim-h0+1):iSim+H-1);
                Q_P_m_H = Q_P_m((iSim-h0+1):iSim+H-1);
                for h = (h0+1):1:(h0+H)
                    [T_z_m_H(h), Q_I_m_H(h), Q_D_m_H(h)] = HVAC(T_z_m_H(h-1), T_z_m_H(h-Model.d_I-1),  T_z_m_H(h-Model.d_D-1), ux, ...
                        T_z_m_H(h-Model.d_I-1), T_z_m_H(h-Model.d_D-1), Q_I_m_H(h-1), Q_D_m_H(h-1), Q_P_m_H(h-1), Model);
                end
                ymH = T_z_m_H(h0+H);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                ePSO = wx - yp;
                em = (yp - T_z_m_H(h0)) + T_z_m_H;
                ymRefH=ym+ePSO*(1-ar^H);

                %                 %             if any(((em - z_sig*sig) < wlb(iSim)) | ((em + z_sig*sig) > wub(iSim)))
                %                 %                 weightYm = 100;
                %                 %             else
                %                 %                 weightYm = 0;
                %                 %             end

                if abs(ux-uZ) > Deltau
                    weightUz =  1e10;
                else
                    weightUz = 0;
                end
                %

                if abs(ux) > Maxu
                    weightU =  1e10;
                else
                    weightU = 0;
                end

                if strcmp(name_test,['HVAC_PFC_test_delta_U_',num2str(Deltau*1e3)])
                    val = (ymH-ymRefH).^2  + lambda_1*(ux.^2) + weightU + weightUz;
                elseif strcmp(name_test,'HVAC_PFC_test')
                    val = (ymH-ymRefH).^2  + lambda_1*(ux.^2) + weightU;
                end

                if val < swarm(i, 4, 1)                 % èe je nova pozicija boljša
                    swarm(i, 3, 1) = swarm(i, 1, 1);    % popravek pozicije x,
                    swarm(i, 3, 2) = swarm(i, 1, 2);    % popravek pozicije x,
                    swarm(i, 4, 1) = val;
                    T_z_m_H_best = T_z_m_H; % najboljša vrednsot delca - pBest
                end
            end

            [temp, gbest] = min(swarm(:, 4, 1));        % gBest
            term = std(swarm(:,1,2)) + std(swarm(:,1,1));


            %-- popravek pozicije delca
            for i = 1 : swarm_size
                swarm(i, 2, :) = rand*inertia*swarm(i, 2, :) + correction_factor*rand*(swarm(i, 3, :) - swarm(i, 1, :)) ...
                    + correction_factor*rand*(swarm(gbest, 3, :) - swarm(i, 1, :));   %x komponenta hitrosti
            end
            %
            %% Plotanje swarma
%             clf
%             plot(swarm(:, 1, 1),swarm(:, 1, 2)-w(iSim), 'x')   %  izris premikanja delca%
%             xlim([-5,5]);
%             ylim([-5,5]);
%              pause(.2)
        end

        if (iter>max_iter)
            disp("Max iter")
        end

        u = swarm(gbest,1,1);
        if (swarm(gbest, 4, 1) >= val_0)
            u = uZ;
        end
        uZ = u;

        wZ = w(iSim); %swarm(gbest,1,2);

        e = wZ - yp;

        %if (wZ + z_sig*sig) < yp
        %   sig = (yp - wZ)/(z_sig);
        %   var2 = (sig^2)*n_var;
        %end

        if iSim == 10
    
            figure(3)
            subplot(2,1,1); hold off;
            plot(0:h0,Yp((end-h0):end),'b'); hold on;
            %plot(h0:h0+H, ypH(h0:h0+H),'b--')
            plot(0:h0,Ym((end-h0):end),'r')
            plot(h0:h0+H,T_z_m_H_best(h0:h0+H),'r--');
            plot(h0:h0+H,ar.^(0:H)*yp+wZ*(1-ar.^(0:H)),'b--');
            plot(0:h0+H,[Ws(iSim-h0:iSim);repmat(wZ,H,1)],'k');
            plot(h0:h0+H,ym+e*(1-ar.^(0:H)),'k--');
            %plot(h0:h0+H,[ym, repelem(ymRefH,1,H)],'m--');
            xline(h0)
            xlabel('Time step')
            ylabel('Temperature [°K]')
            xlim('tight')
            xticks(k_tics);
            xticklabels(labels);
            legend('Process past','Model past','Model future','Reference value',...
                'Reference model','Desired model trajectory','location','best',...
                'FontSize',8)
            pause(0.0001)

            name = ['..\HVAC_Images\',name_test,'.pdf'];
            exportgraphics(gcf,name,'BackgroundColor','none');
            pause(0.0001)


        end

        % simulacija procesa
        d = -0.2; noise = 0; M.N;
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
            T_z_p(iSim-M.d_I-1), T_z_p(iSim-M.d_D-1), Q_I_p(iSim-1), Q_D_p(iSim-1), Q_P_p(iSim-1), M);
        yp = T_z_p(iSim) + noise*randn(1,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if abs(w(iSim)- w(iSim-1))>0
            %w_sigma = w(iSim);
            var2 = sigma_0;
            n_var = n_var_0;
            mu = yp;
            T_z_m((end-h0):end) = T_z_p((end-h0):end);

        end

        e_var = yp - mu; %Distance of new data from center
        mu = mu + 1/(1 + n_var)*e_var; %Center update
        var2  = var2 + e_var*(yp - mu)'; %un-normalized covariance matrix
        n_var = n_var + 1; %Increase number of samples in cluster
        sig = sqrt(var2/n_var);

        z_sig = tinv(0.95, n_var);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % simulacije modela
        um = u;
        % ym=am*ymZ+bm*um;
        % ymZ=ym;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [T_z_m(iSim), Q_I_m(iSim), Q_D_m(iSim)] = HVAC(T_z_m(iSim-1), T_z_m(iSim-Model.d_I-1), T_z_m(iSim-Model.d_D-1), um, ...
            T_z_m(iSim-Model.d_I-1), T_z_m(iSim-Model.d_D-1), Q_I_m(iSim-1), Q_D_m(iSim-1), Q_P_m(iSim-1), Model);
        ym = T_z_m(iSim);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Yp = [Yp; yp];
        Ym = [Ym; ym];
        U = [U; up];
        W = [W; w(iSim)];
        Ws = [Ws; wZ];
        V = [V; sig];
        Z = [Z;z_sig];
        if (mod(iSim,100) == 0)

            h1 = figure(1);
            subplot(3,1,1:2); hold off;
            p11 = plot(Yp,'b'); hold on;
            p12 = plot(Yp + z_sig*V,'-.',"color",color(1,:));
            plot(Yp - z_sig*V,'-.',"color",color(1,:));
            p3 = plot(Ym,'r');
            p4 = plot(wub,'k');
            plot(wlb,'k')
            p5 = plot(w,'k--');

            legend([p11,p12,p3,p4,p5],'System output','Confidence interval',...
                'Model output','Reference interval','Reference value','location','best')
            xlim('tight')
            xlabel('Time step')


            h2 = figure(2);
            subplot(3,1,1); hold off;
            h2 = plot(U,'b');
            ylabel({'Control','signal'})
            xlim('tight')

            subplot(3,1,2); hold off;
            plot((Q_D_p(1:length(U)-1) + M.act_P*Q_P_p(1:length(U)-1)'),"color",color(3,:)); hold on;
            plot((Q_I_p(1:length(U)-1)),'r');
            legend('Disturbance', 'Control','location','best')
            ylabel({'Thermal','energy'})
            xlim('tight')
            subplot(3,1,3); hold off;
            plot(u_sum,"color",color(1,:));
            ylabel('Cost')
            legend('Sum of squared inputs','Location','best')
            xlabel('Time step')
            xlim('tight')
            pause(0.001)


        end
    end
    t2 = clock;
    save(name_test)

end

load HVAC_PFC_test

Povprecni_Cas_Enega_Izracuna = etime(t2,t1)/Nsim;

h1 = figure(4);
subplot(3,1,1:2); hold off;
p11 = plot(w,'k--'); hold on;
p12 = plot(Yp,'b'); 
%p3 = plot(Ym,'b--');
ylim('tight')
ylabel('Temperature')
subplot(3,1,3); hold off;
p21 = plot(U,'b');hold on;
ylabel({'Control','signal'})
xlim('tight')
xlabel('Time step')

load HVAC_PFC_test_delta_U_10
h1 = figure(4);
subplot(3,1,1:2);
p14 = plot(Yp,'r'); 
ylim('tight')
ylabel('Temperature')
subplot(3,1,3);
p22 = plot(U,'r');
ylabel({'Control','signal'})
xlim('tight')
xlabel('Time step')

load HVAC_PFC_test_delta_U_1
h1 = figure(4);
subplot(3,1,1:2);
p15 = plot(Yp,'m'); hold on;
ylim('tight')
ylabel('Temperature')
legend('Reference value','Control law without ${\Delta}u$','Control law with ${\Delta}u = 0.01$',...
    'Control law with ${\Delta}u = 0.001$','location', 'best','interpreter','latex',...
                'FontSize',10)
xlim('tight') 
subplot(3,1,3); hold on;
p23 = plot(U,'m');
ylabel({'Control','signal'})
xlim('tight')
xlabel('Time step')
% legend('Control law without ${\Delta}u$','Control law with ${\Delta}u = 0.01$',...
%     'Control law with ${\Delta}u = 0.001$','location', 'best','interpreter','latex',...
%                 'FontSize',10)

figure(4);
name = '..\HVAC_Images\HVAC_PFC_test.pdf';
exportgraphics(gcf,name,'BackgroundColor','none');
pause(0.0001)

% figure(2);
% name = '..\HVAC_Images\HVAC_PFC_outcov_cost.pdf';
% exportgraphics(gcf,name,'BackgroundColor','none');
% pause(0.0001)
