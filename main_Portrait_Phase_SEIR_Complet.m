function main_Portrait_Phase_SEIR_Complet()
    % =====================================================================
    % PARAMÈTRES BIOLOGIQUES ET DÉMOGRAPHIQUES FIXES
    % =====================================================================
    Lambda = 10;     % Taux de recrutement
    alpha  = 0.5;    % Taux d'incubation (E -> I)
    gamma  = 0.3;    % Taux de guérison
    mu     = 0.1;    % Taux de mortalité naturelle
    
    % Temps de simulation pour les trajectoires
    tspan = [0 200]; 
    
    % Conditions initiales diversifiées pour bien remplir l'espace de phase (S, I)
    % Format : [S0, E0, I0, R0]
    S0_samples = [20,  50, 150, 200,  20, 80, 180, 240];
    I0_samples = [80,  90,  80,  60,  10,  5,  15,  5];
    E0_default = 10;
    R0_default = 0;
    
    % Configuration de la fenêtre d'affichage
    figure('Position', [50, 100, 1350, 600], 'Color', [1 1 1]);
    
    % =====================================================================
    % SOUS-GRAPHIQUE 1 : CAS R0 < 1 (beta = 0.0002) -> EXTINCTION
    % =====================================================================
    beta_cas1 = 0.0002; 
    
    % Calcul du R0 pour le cas 1
    R0_1 = (beta_cas1 * Lambda * alpha) / (mu * (mu + alpha) * (mu + gamma));
    
    % Point d'équilibre sans maladie (DFE)
    S_DFE1 = Lambda / mu; % Égal à 100
    I_DFE1 = 0;
    
    subplot(1, 2, 1);
    hold on; box on; grid on;
    
    % 1. Champ de directions (Vecteurs normalisés pour le sens du mouvement)
    [S_mesh, I_mesh] = meshgrid(linspace(5, 250, 20), linspace(1, 100, 20));
    dS1 = Lambda - beta_cas1.*S_mesh.*I_mesh - mu.*S_mesh;
    dI1 = (alpha/(mu+alpha)).*beta_cas1.*S_mesh.*I_mesh - (mu+gamma).*I_mesh; 
    Norm1 = sqrt(dS1.^2 + dI1.^2);
    quiver(S_mesh, I_mesh, dS1./Norm1, dI1./Norm1, 0.4, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    
    % 2. Tracé des trajectoires (Flux dynamique)
    for k = 1:length(S0_samples)
        X0 = [S0_samples(k); E0_default; I0_samples(k); R0_default];
        [~, X] = ode45(@(t,X) seir_system(t,X,Lambda,beta_cas1,alpha,gamma,mu), tspan, X0);
        
        % Trajectoire en bleu (Extinction)
        plot(X(:,1), X(:,3), 'Color', [0 0.4 0.8], 'LineWidth', 2);
        % Point de départ de la trajectoire
        plot(X(1,1), X(1,3), 'ok', 'MarkerFaceColor', [0 0.4 0.8], 'MarkerSize', 5);
    end
    
    % 3. Tracé du point d'équilibre DFE (STABLE GLOBALE)
    plot(S_DFE1, I_DFE1, 'Marker', 'h', 'MarkerSize', 14, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0.8 0]);
    
    % Mises en page du graphique 1
    xlabel('Individus Sains (S)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Individus Infectés (I)', 'FontSize', 11, 'FontWeight', 'bold');
    title({['Portrait de Phase pour R_0 = ', numless(R0_1), ' \le 1'], ...
           'Comportement : Toutes les trajectoires finissent au point DFE (Extinction)'}, 'FontSize', 11);
    xlim([0 260]); ylim([0 110]);
    legend('Sens du mouvement (Champ)', 'Trajectoires', 'Départ (t=0)', 'DFE Stable Global (S^*, 0)', 'Location', 'northeast');
    
    % =====================================================================
    % SOUS-GRAPHIQUE 2 : CAS R0 > 1 (beta = 0.008) -> PERSISTANCE ENDÉMIQUE
    % =====================================================================
    beta_cas2 = 0.005; 
    
    % Calcul du R0 pour le cas 2
    R0_2 = (beta_cas2 * Lambda * alpha) / (mu * (mu + alpha) * (mu + gamma));
    
    % Point d'équilibre sans maladie (DFE) -> Devenu INSTABLE
    S_DFE2 = Lambda / mu; 
    I_DFE2 = 0;
    
    % Point d'équilibre endémique (EE) -> STABLE GLOBALE
    S_EE = ((mu + alpha) * (mu + gamma)) / (alpha * beta_cas2); % Égal à 48
    I_EE = (Lambda - mu * S_EE) / (beta_cas2 * S_EE);          % Égal à 21.66
    
    subplot(1, 2, 2);
    hold on; box on; grid on;
    
    % 1. Champ de directions
    dS2 = Lambda - beta_cas2.*S_mesh.*I_mesh - mu.*S_mesh;
    dI2 = (alpha/(mu+alpha)).*beta_cas2.*S_mesh.*I_mesh - (mu+gamma).*I_mesh;
    Norm2 = sqrt(dS2.^2 + dI2.^2);
    quiver(S_mesh, I_mesh, dS2./Norm2, dI2./Norm2, 0.4, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
    
    % 2. Tracé des trajectoires
    for k = 1:length(S0_samples)
        X0 = [S0_samples(k); E0_default; I0_samples(k); R0_default];
        [~, X] = ode45(@(t,X) seir_system(t,X,Lambda,beta_cas2,alpha,gamma,mu), tspan, X0);
        
        % Trajectoire en rouge/orange (Persistance endémique)
        plot(X(:,1), X(:,3), 'Color', [0.8 0.2 0], 'LineWidth', 2);
        % Point de départ
        plot(X(1,1), X(1,3), 'ok', 'MarkerFaceColor', [0.8 0.2 0], 'MarkerSize', 5);
    end
    
    % 3. Tracé des deux points d'équilibre
    % Le DFE qui est devenu INSTABLE (marqué en ROUGE)
    plot(S_DFE2, I_DFE2, 'Marker', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0 0]); 
    % L'équilibre Endémique (EE) qui est STABLE GLOBAL (marqué en ÉTOILE VERTE)
    plot(S_EE, I_EE, 'Marker', 'pentagram', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0 0.8 0]);
    
    % Mises en page du graphique 2
    xlabel('Individus Sains (S)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Individus Infectés (I)', 'FontSize', 11, 'FontWeight', 'bold');
    title({['Portrait de Phase pour R_0 = ', numless(R0_2), ' > 1'], ...
           'Comportement : Le DFE est fui (instable), convergence vers l''Équilibre Endémique'}, 'FontSize', 11);
    xlim([0 260]); ylim([0 110]);
    legend('Sens du mouvement (Champ)', 'Trajectoires', 'Départ (t=0)', 'DFE Instable (S^*, 0)', 'EE Stable Global (S^*, I^*)', 'Location', 'northeast');
end

% =====================================================================
% ENCAPSULATION DES ÉQUATIONS DIFFÉRENTIELLES DU SYSTÈME SEIR
% =====================================================================
function dXdt = seir_system(~, X, Lambda, beta, alpha, gamma, mu)
    S = X(1); E = X(2); I = X(3); R = X(4);
    
    dSdt = Lambda - beta*S*I - mu*S;
    dEdt = beta*S*I - (mu + alpha)*E;
    dIdt = alpha*E - (mu + gamma)*I;
    dRdt = gamma*I - mu*R;
    
    dXdt = [dSdt; dEdt; dIdt; dRdt];
end

% Fonction auxiliaire pour formater l'affichage de R0 dans les titres
function str = numless(val)
    str = num2str(val, '%.2f');
end