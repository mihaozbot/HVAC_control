function  [Yp, Ym, U, W, Ws, V, Z, u_sum, Q_D_p, Q_I_p, wlb, wub] = compute_HVAC_PFC_interval(Deltau, lambda_1, lambda_2, ar, H, alpha, n_var_0, var_0)

color = getcolors(7,"default");

HVAC_par %Load model
Q_P = Q_P_data;

c_P = 1; %People inforamtion active
c_O = 1; %Outside inforamtion active
k_0 = max([M.d_I,M.d_D])+2;

w = 295*ones(10000,1);
w(2000:5000) = 296;
w(8000:10000) = 294;
w = w(1:1:end);
Nsim = length(w);
dw = 2;
wlb = w - dw;
wub = w + dw;

w_shift = 200;
wlb_cor = min([[wlb(1:w_shift);wlb(1:end-w_shift)],...
    [wlb(w_shift+1:end);wlb((end-w_shift+1):end)]],[],2);
wub_cor = max([[wub(1:w_shift);wub(1:end-w_shift)],...
    [wub(w_shift+1:end);wub((end-w_shift+1):end)]],[],2);

%  proces 1. red
Q_I_p = 0*ones(k_0-1,1);
Q_D_p = 0*ones(k_0-1,1);
T_z_p = w(k_0-1)*ones(k_0-1,1);
T_O_p = c_O*T_o_data + (1-c_O)*295;
Q_P_p  = Q_P;

%  model procesa 1. reda
Model = M;
Model.N = 0; %0.5000;
Model.K_D = 0.7*M.K_D;  %11.0657;
Model.tau_D = 0.7*M.tau_D;  %1200;
Model.b_I = 0.7*M.b_I;  %-117.2899;

Q_I_m = 0*ones(k_0-1,1);
Q_D_m = 0*ones(k_0-1,1);
T_z_m = w(k_0-1)*ones(k_0-1,1);
T_O_m = c_O*T_o_data + (1-c_O)*295;
Q_P_m = c_P*Q_P;

M.act_P = 1;
Model.act_P = 0; %0;

% parametri simulacije
yp = T_z_p(end); ym = T_z_m(end);

n_var = n_var_0;
z_sig =  tinv(alpha, n_var);
var = var_0;
sig = sqrt(var/n_var);
mu = yp; %(wlb(1) + 3*sqrt(sigma2/n_var));

U = 0*ones(k_0-1,1); Ym = T_z_m; Yp = T_z_p;  W = T_z_p;  uZ=0; Ws = T_z_p;
V = sig*ones(k_0-1,1);
Z = 2*ones(k_0-1,1);

% referenèni model
br = (1-ar);
N1 = 1;
N2 = H;
u_sum = 0*ones(k_0-1,1);

% ---- inicializacija -----
swarm_size = 10;
for i = 1 : swarm_size
    swarm(i, 1, 1) = i;
end

u = 0; Maxu = 4;
t1 = clock;

h0 = max([M.d_I,M.d_D])+2-1;
k_tics = 0:floor((h0+H)/10):(h0+H);

T_z_m_H = zeros(h0+H,1);
Q_I_m_H = zeros(h0+H,1);
Q_D_m_H = zeros(h0+H,1);

for iSim = k_0:Nsim

    % Insert the break condition after Yp, Z, and V are updated
    % try
    %     if (any((Yp + Z .* V ) > wub_cor(1:length(Yp))) + ...
    %         any((Yp - Z .* V) < wlb_cor(1:length(Yp))) + ...
    %         ((2*best_cost)<u_sum(end))) > 0
    %         u_sum = inf;
    %         return;  % Exit the function early
    %     end
    % catch
    %     any((Yp + Z .* V ) > wub_cor(1:length(Yp)))
    %     disp("Error")
    %     disp(size((Yp+Z.* V)))
    %     disp(size(wub_cor(1:length(Yp))))
    %     disp(((2*best_cost)<u_sum(end)))
    % 
    % end

    inertia = 1.0;
    correction_factor = 2;

    % ---- inicializacija -----
    for i = 1 : swarm_size
        %      swarm(i, 1, 1) =  u - round(swarm_size/2)*Maxu/swarm_size + i * (Deltau);
        swarm(i, 1, 1) = Deltau/2- i * (Deltau/swarm_size) + u ;
    end

    swarm(:, 1, 2) =  repmat(floor( 1:swarm_size-swarm_size/2).*((wub(iSim)-wlb(iSim))/swarm_size)*1 + w(iSim),1,2);
    %         swarm((swarm(:, 1, 2) < min(wlb(iSim) + z_sig*sig, w(iSim))) |...
    %             (swarm(:, 1, 2)  > max(wub(iSim) - z_sig*sig,w(iSim))), 1, 2) = w(iSim);


    swarm(:,3,1) = swarm(:,1,1);
    swarm(:,3,2) = swarm(:,1,2);
    val_0 = 1e5;
    swarm(:, 4, 1) = 1e5;          % najboljša vrednsot do sedaj
    swarm(:, 4, 2) = NaN;          % najboljša vrednsot do sedaj
    swarm(:, 2, :) = 0;             % zaèetna hitrost

    term = 1000; termMax = 1e-4;

    % izdelaj termination criteria
    max_iter = 100;
    iter = 0;
    eb_error = 0;
    error_weight = 1e2;
    while ((term > termMax) && (iter<max_iter))

        %for iter = 1 : iterations
        iter = iter +1;

        %-- evaluacija pozicije in kriterijske funkcije ---
        for i = 1 : swarm_size

            swarm(i, 1, :) = swarm(i, 1, :) + swarm(i, 2,:); % popravek pozicije x
            ux = swarm(i, 1, 1);
            wx = swarm(i, 1, 2);

            if abs(ux-uZ) > Deltau
                weightUz =  error_weight;
            else
                weightUz = 0;
            end
            %

            if abs(ux) > Maxu
                weightU =  error_weight;
            else
                weightU = 0;
            end

            if (weightU + weightUz) > 0
                continue
            end

            if lambda_1*(ux.^2) > swarm(i, 4, 1)
                continue
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            T_z_m_H(1:h0) = T_z_m((end-h0+1):end);
            Q_I_m_H(1:h0) = Q_I_m((end-h0+1):end);
            Q_D_m_H(1:h0) = Q_D_m((end-h0+1):end);
            T_O_m_H = T_O_m((iSim-h0+1):iSim+H-1);
            Q_P_m_H = Q_P_m((iSim-h0+1):iSim+H-1);
            for h = (h0+1):1:(h0+H)
                [T_z_m_H(h), Q_I_m_H(h), Q_D_m_H(h)] = HVAC(T_z_m_H(h-1), T_z_m_H(h-Model.d_I-1),  T_z_m_H(h-Model.d_D-1), ux, ...
                    T_O_m_H(h-Model.d_I-1), T_O_m_H(h-Model.d_D-1), Q_I_m_H(h-1), Q_D_m_H(h-1), Q_P_m_H(h-1), Model);
            end
            ymH = T_z_m_H(h0+H);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ePSO = wx - yp;
            ypH =  yp + (T_z_m_H - T_z_m_H(h0));
            ymRefH = ym+ePSO*(1-ar^H);


            elb = w(iSim)-(ypH(end) - z_sig*sig);
            elb = max(elb, 0);
            eub = (ypH(end) + z_sig*sig)-w(iSim);
            eub = max(eub, 0);
            eb = elb^2 + eub^2;
            
            if ( any((Yp + Z .* V ) > wub_cor(1:length(Yp))) + ...
             any((Yp - Z .* V) < wlb_cor(1:length(Yp))) )>0
                eb = eb + error_weight;
            end

            val = (ymH-ymRefH).^2  + lambda_1*(ux.^2) + lambda_2*sum(eb.^2); % izraèun kriterijske funkcije

            if val < swarm(i, 4, 1)                 % če je nova pozicija boljša
                swarm(i, 3, 1) = swarm(i, 1, 1);    % popravek pozicije x,
                swarm(i, 3, 2) = swarm(i, 1, 2);    % popravek pozicije x,
                swarm(i, 4, 1) = val;

            end
        end

        [~, gbest] = min(swarm(:, 4, 1)); % gBest

        term = std(swarm(:,1,2)) + std(swarm(:,1,1));

        %-- popravek pozicije delca
        for i = 1 : swarm_size
            swarm(i, 2, :) = rand*inertia*swarm(i, 2, :) + correction_factor*rand*(swarm(i, 3, :) - swarm(i, 1, :)) ...
                + correction_factor*rand*(swarm(gbest, 3, :) - swarm(i, 1, :));   %x komponenta hitrosti
        end

    end

    if (iter>max_iter)
        disp("Max iter")
    end

    u = swarm(gbest,1,1);
    %         if (swarm(gbest, 4, 1) >= val_0)
    %             u = uZ;
    %         end
    uZ = u;
    wZ = swarm(gbest,1,2);
    e = wZ - yp;

    % simulacija procesa
    d = -0.2; noise = M.N;
    up = u;
    u_sum(iSim) = u_sum(iSim-1)+up^2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [T_z_p(iSim), Q_I_p(iSim), Q_D_p(iSim)] = HVAC(T_z_p(iSim-1), T_z_p(iSim-M.d_I-1), T_z_p(iSim-M.d_D-1), up,...
        T_O_p(iSim-M.d_I-1), T_O_p(iSim-M.d_D-1), Q_I_p(iSim-1), Q_D_p(iSim-1), Q_P_p(iSim-1), M);
    yp = T_z_p(iSim) + noise*randn(1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if abs(w(iSim)- w(iSim-1))>0
        %w_sigma = w(iSim);
        var = var_0;
        n_var = n_var_0;
        mu = yp;
        T_z_m((end-h0):end) = T_z_p((end-h0):end);

    end

    e_var = yp - mu; %Distance of new data from center
    mu = mu + 1/(1 + n_var)*e_var; %Center update
    var  = var + e_var*(yp - mu)'; %un-normalized covariance matrix
    n_var = n_var + 1; %Increase number of samples in cluster
    sig = sqrt(var/n_var);

    z_sig = tinv(alpha, n_var);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % simulacije modela
    um = u;

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
    Z = [Z; z_sig];

    if (mod(iSim,10000) == 0)

        h1 = figure(1);
        set(h1,'Position',[680,558,2*560,420])
        subplot(3,1,1:2); hold off;
        p1 = plot(Yp,'b'); hold on;
        p2 = plot(Yp + Z.*V,'-.',"color",color(1,:));
        plot(Yp - Z.*V,'-.',"color",color(1,:));
        p3 = plot(Ym,'r');
        p4 = plot(wub,'k');
        plot(wlb,'k')
        legend([p1,p2,p3,p4],'System output','Confidence interval','Model output','Reference interval','location','best')
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
        pause(0.00001)

    end


end


end
