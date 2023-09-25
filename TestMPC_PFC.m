% MPC realiziran z PSO
% 1.3.2013

clear all
close all


%  proces 1. red
ap=0.95; bp=0.03; uZ=0;

%  model procesa 1. reda
am = 1 * ap; bm = 0.9 * bp; Km=bm/(1-am);

% parametri simulacije
Nsim = 800;
ypZ=0; yp=0; ymZ=0; ym=0;
w=0.75;
U=[]; Ym=[]; Yp=[]; Yp2=[]; W=[];

% referenèni model
ar=0.7; br=(1-ar); H=5;


% ---- inicializacija -----
swarm_size = 5;
for i = 1 : swarm_size
    swarm(i, 1, 1) = i;
end

u = 0; Maxu = 1.5; Deltau = 0.075;

t1 = clock;
for iSim=1:Nsim
    
    if iSim > 400 
        w = -1;
    end
    
    e=w-yp;
    ymRefH=ym+e*(1-ar^H);    

    % PSO algoritem
    iterations = 100;
    inertia = 1.0;
    correction_factor = 2;

    % ---- inicializacija -----
    for i = 1 : swarm_size
 %      swarm(i, 1, 1) =  u - round(swarm_size/2)*Maxu/swarm_size + i * (Deltau); 
       swarm(i, 1, 1) =  i * (Maxu/swarm_size)*0.1 + u ; 
    end

    swarm(:, 4, 1) = 1000;          % najboljša vrednsot do sedaj
    swarm(:, 2, :) = 0;             % zaèetna hitrost
    
    term = 1000; termMax = 0.02;
    % izdelaj termination criteria
    while term > termMax
    %for iter = 1 : iterations

        %-- evaluacija pozicije in kriterijske funkcije ---
        for i = 1 : swarm_size
            swarm(i, 1, 1) = swarm(i, 1, 1) + swarm(i, 2, 1);     % popravek pozicije x 
            ux = swarm(i, 1, 1);
       
            Phi = [ym ux]; Th = [am bm]';
            ymH = ModelBasedPrediction_MeanLevel(Phi,Th,H);
 
            if ux > Maxu, 
                weightU = 1000; 
            else
                weightU = 0;
            end

            if abs(ux-uZ) > Deltau,
                weightUz = 1000; 
            else
                weightUz = 0;
            end

            val= (ymH-ymRefH).^2 + weightU + weightUz;        % izraèun kriterijske funkcije

           
            
            if val < swarm(i, 4, 1)                 % èe je nova pozicija boljša
                swarm(i, 3, 1) = swarm(i, 1, 1);    % popravek pozicije x,
                swarm(i, 4, 1) = val;               % najboljša vrednsot delca - pBest
            end
        end
        
        swarm(1:swarm_size,1,1);
        term = std(swarm(1:swarm_size,1,1));
    
        [temp, gbest] = min(swarm(:, 4, 1));        % gBest

        %-- popravek pozicije delca
        for i = 1 : swarm_size
            swarm(i, 2, 1) = rand*inertia*swarm(i, 2, 1) + correction_factor*rand*(swarm(i, 3, 1) - swarm(i, 1, 1)) + correction_factor*rand*(swarm(gbest, 3, 1) - swarm(i, 1, 1));   %x komponenta hitrosti
        end

        %% Plotanje swarma
    %    clf
    %    plot(swarm(:, 1, 1), 'x')   %  izris premikanja delca
    %    axis([0 10 -3 3])
    %    pause(.2)
    end

    
    u=swarm(gbest,1,1);
    uZ=u;

    % simulacija procesa
    d=-0.2; noise=0.005;
    up=u;
    if iSim > round(2*Nsim/4)
        up=u+d;
    end
    
%     yp = ap*ypZ+bp*up;
%     yp=yp+noise*randn(1,1);
%     ypZ=yp;
    
    Phi_p = [ypZ up]; Th_p = [ap bp]';
    ypZ = ModelBasedPrediction_MeanLevel(Phi_p,Th_p,1);
    yp = ypZ + noise*randn(1,1);
  
    % simulacije modela
    um=u;
%     ym=am*ymZ+bm*um;
%     ymZ=ym;

    Phi_m = [ym um]; Th_m = [am bm]';
    ym = ModelBasedPrediction_MeanLevel(Phi_m,Th_m,1);
   
    
    Yp=[Yp; yp];
    Ym=[Ym; ym];
    U=[U; up];
    W=[W; w];
    
end
t2 = clock;


Povprecni_Cas_Enega_Izracuna = etime(t2,t1)/Nsim

plot([W Yp Ym U])
legend('reference','process','model','input')
    
    
