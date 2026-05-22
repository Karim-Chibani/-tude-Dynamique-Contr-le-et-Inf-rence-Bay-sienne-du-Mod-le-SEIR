clc; clear; close all;

%% =====================================================================
%% 1. PARAMÈTRES ET STRUCTURE GLOBALE
%% =====================================================================
% On regroupe les paramètres dans une structure globale 'p' pour éviter 
% les redéfinitions locales dans les fonctions (Source d'erreurs).
p.Lambda = 10;
p.beta   = 0.008;
p.alpha  = 0.5;
p.gamma  = 0.3;
p.mu     = 0.1;
p.rho    = 1;     % Poids du coût du contrôle

% Conditions initiales des variables d'état
S0 = 500; E0 = 10; I0 = 5; R0 = 0;
X0 = [S0; E0; I0; R0];

% Paramètres temporels et grille
T = 200;
tspan = [0 T];
N = 300;                  % Augmenté pour une meilleure précision d'intégration
tgrid = linspace(0, T, N);

% Initialisation du contrôle (première estimation : aucun contrôle)
u = zeros(1, N); 

% Paramètres de l'algorithme itératif
max_iter = 100;
tol = 1e-4;

%% =====================================================================
%% 2. ALGORITHME FORWARD-BACKWARD SWEEP
%% =====================================================================
fprintf('Début des itérations du contrôle optimal...\n');

for iter = 1:max_iter
    % Sauvegarde de l'ancien contrôle pour tester la convergence
    u_old = u;
    
    % --- ÉTAPE A : ÉQUATION DIRECTE (FORWARD SYSTEM) ---
    % On résout le système SEIR du passé vers le futur
    [t_fwd, X_fwd] = ode45(@(t,X) SEIR_forward(t, X, tgrid, u, p), tspan, X0);
    
    % Interpolation des états sur la grille fixe 'tgrid' pour le système adjoint
    S_interp = interp1(t_fwd, X_fwd(:,1), tgrid, 'linear', 'extrap');
    E_interp = interp1(t_fwd, X_fwd(:,2), tgrid, 'linear', 'extrap');
    I_interp = interp1(t_fwd, X_fwd(:,3), tgrid, 'linear', 'extrap');
    X_interp = [S_interp', E_interp', I_interp']; % Matrice Nx3 des états interpolés

    % --- ÉTAPE B : ÉQUATION ADJOINTE (BACKWARD SYSTEM) ---
    % Condition finale transversale : lambda(T) = [0, 0, 0, 0]
    lambdaT = [0; 0; 0; 0];
    
    % Résolution du futur vers le passé (flip de tspan)
    [t_bwd, lambda_bwd] = ode45(@(t,lam) adjoint_backward(t, lam, tgrid, X_interp, u, p), flip(tspan), lambdaT);
    
    % Re-synchronisation temporelle des multiplicateurs de Lagrange (du présent vers le futur)
    lambda_interp = zeros(N, 4);
    for k = 1:4
        lambda_interp(:, k) = interp1(t_bwd, lambda_bwd(:, k), tgrid, 'linear', 'extrap');
    end

    % --- ÉTAPE C : MISE À JOUR ET PROJECTION DU CONTRÔLE ---
    u_new = zeros(1, N);
    for i = 1:N
        % Condition d'optimalité dH/du = 0
        u_new(i) = (I_interp(i) * (lambda_interp(i, 3) - lambda_interp(i, 4))) / p.rho;
        
        % Projection stricte sur les bornes admissibles [0, 1]
        u_new(i) = max(0, min(1, u_new(i)));
    end
    
    % --- ÉTAPE D : RELAXATION ET CRITÈRE D'ARRÊT ---
    % Combinaison convexe (0.5/0.5) pour éviter les oscillations numériques numériques
    u = 0.5 * u + 0.5 * u_new;
    
    % Vérification du critère de convergence
    erreur = sum(abs(u - u_old)) / sum(u);
    if (erreur < tol) && (iter > 2)
        fprintf('Convergence atteinte à l''itération %d ! (Erreur = %e)\n', iter, erreur);
        break;
    end
end

%% =====================================================================
%% 3. AFFICHAGE DES RÉSULTATS (VISUALISATION RIGOUREUSE)
%% =====================================================================
figure('Position', [100, 100, 1100, 450]);

% Subplot 1 : Trajectoire des infectés interpolés
subplot(1, 2, 1);
plot(tgrid, I_interp, 'r-', 'LineWidth', 2.5);
xlabel('Temps (jours)', 'FontSize', 11);
ylabel('Population Infectée I(t)', 'FontSize', 11);
title('Dynamique des Infectés sous Contrôle Optimal', 'FontSize', 12);
grid on;

% Subplot 2 : Profil du contrôle optimal obtenu
subplot(1, 2, 2);
plot(tgrid, u, 'b-', 'LineWidth', 2.5);
xlabel('Temps (jours)', 'FontSize', 11);
ylabel('Effort du contrôle u(t)', 'FontSize', 11);
title('Profil de la Stratégie de Contrôle Optimale', 'FontSize', 12);
ylim([-0.05 1.05]);
grid on;

%% =====================================================================
%% 4. FONCTIONS NUMÉRIQUES INTERNES
%% =====================================================================

function dX = SEIR_forward(t, X, tgrid, u_vec, p)
    S = X(1); E = X(2); I = X(3); R = X(4);
    
    % Interpolation de la valeur instantanée du contrôle à l'instant t
    u_t = interp1(tgrid, u_vec, t, 'linear', 'extrap');
    u_t = max(0, min(1, u_t)); % Sécurité des bornes
    
    % Système différentiel direct
    dS = p.Lambda - p.beta * S * I - p.mu * S;
    dE = p.beta * S * I - (p.mu + p.alpha) * E;
    dI = p.alpha * E - (p.mu + p.gamma) * I - u_t * I;
    dR = p.gamma * I + u_t * I - p.mu * R;
    
    dX = [dS; dE; dI; dR];
end

function dlambda = adjoint_backward(t, lambda, tgrid, X_interp, u_vec, p)
    lambdaS = lambda(1);
    lambdaE = lambda(2);
    lambdaI = lambda(3);
    lambdaR = lambda(4);
    
    % Récupération de l'état interpolé instantané au temps t courant
    S_t = interp1(tgrid, X_interp(:, 1), t, 'linear', 'extrap');
    I_t = interp1(tgrid, X_interp(:, 3), t, 'linear', 'extrap');
    
    % Récupération du contrôle instantané au temps t courant
    u_t = interp1(tgrid, u_vec, t, 'linear', 'extrap');
    u_t = max(0, min(1, u_t));
    
    % Système adjoint linéaire rétrograde (issu de la matrice jacobienne du Hamiltonien)
    dlambdaS = lambdaS * (p.beta * I_t + p.mu) - lambdaE * (p.beta * I_t);
    
    dlambdaE = lambdaE * (p.mu + p.alpha) - lambdaI * p.alpha;
    
    % Le terme "-1" provient de la minimisation de l'intégrale du nombre d'infectés dans le coût J
    dlambdaI = -1 + lambdaS * (p.beta * S_t) - lambdaE * (p.beta * S_t) ...
               + lambdaI * (p.mu + p.gamma + u_t) - lambdaR * (p.gamma + u_t);
               
    dlambdaR = p.mu * lambdaR;
    
    dlambda = [dlambdaS; dlambdaE; dlambdaI; dlambdaR];
end