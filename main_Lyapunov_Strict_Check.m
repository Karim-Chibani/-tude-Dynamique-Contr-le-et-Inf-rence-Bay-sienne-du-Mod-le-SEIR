function main_Lyapunov_Strict_Check()
    % --- Paramètres démographiques et biologiques fixes ---
    Lambda = 10;
    alpha  = 0.5;
    gamma  = 0.3;
    mu     = 0.1;
    
    tspan = [0 150]; 
    X0 = [200; 30; 50; 10]; % Conditions initiales [S0, E0, I0, R0]
    
    % Configuration de la figure
    figure('Position', [100, 100, 1200, 600], 'Color', [1 1 1]);
    
    % =====================================================================
    % CAS 1 : R0 <= 1 (Extinction) -> beta = 0.0002
    % =====================================================================
    beta1 = 0.0002; 
    S_star1 = Lambda / mu; % Point d'équilibre S_DFE
    
    [t1, X1] = ode45(@(t,X) seir_eq(t,X,Lambda,beta1,alpha,gamma,mu), tspan, X0);
    [V1, dV1] = calculate_lyapunov(X1, Lambda, beta1, alpha, gamma, mu, S_star1);
    
    % Tracé Cas 1 : V(t) et dV/dt
    subplot(2, 2, 1);
    plot(t1, V1, 'g-', 'LineWidth', 2.5); grid on;
    title('V(t) pour R_0 \le 1 (Strictement Décroissante)', 'FontSize', 10);
    ylabel('Valeur de V(S,E,I)');
    
    subplot(2, 2, 3);
    plot(t1, dV1, 'r-', 'LineWidth', 2); grid on;
    line([0 tspan(2)], [0 0], 'Color', 'k', 'LineStyle', '--'); % Ligne 0
    title('dV/dt pour R_0 \le 1 (Toujours Négative \le 0)', 'FontSize', 10);
    xlabel('Temps (jours)'); ylabel('Axe dV/dt');
    
    % =====================================================================
    % CAS 2 : R0 > 1 (Persistance Endémique) -> beta = 0.005
    % =====================================================================
    beta2 = 0.005;
    % Pour R0 > 1, S* devient le S de l'équilibre Endémique (S_EE)
    S_star2 = ((mu + alpha)*(mu + gamma)) / (alpha * beta2); 
    
    [t2, X2] = ode45(@(t,X) seir_eq(t,X,Lambda,beta2,alpha,gamma,mu), tspan, X0);
    [V2, dV2] = calculate_lyapunov(X2, Lambda, beta2, alpha, gamma, mu, S_star2);
    
    % Tracé Cas 2 : V(t) et dV/dt
    subplot(2, 2, 2);
    plot(t2, V2, 'b-', 'LineWidth', 2.5); grid on;
    title('V(t) pour R_0 > 1 (Strictement Décroissante)', 'FontSize', 10);
    ylabel('Valeur de V(S,E,I)');
    
    subplot(2, 2, 4);
    plot(t2, dV2, 'r-', 'LineWidth', 2); grid on;
    line([0 tspan(2)], [0 0], 'Color', 'k', 'LineStyle', '--');
    title('dV/dt pour R_0 > 1 (Toujours Négative \le 0)', 'FontSize', 10);
    xlabel('Temps (jours)'); ylabel('Axe dV/dt');
end

% --- Fonction pour calculer V et sa dérivée dV/dt à chaque instant ---
function [V, dVdt] = calculate_lyapunov(X, Lambda, beta, alpha, gamma, mu, S_star)
    S = X(:, 1);
    E = X(:, 2);
    I = X(:, 3);
    
    % 1. Calcul de la fonction V(t)
    V = (S - S_star - S_star .* log(S ./ S_star)) + E + ((mu + alpha) / alpha) .* I;
    
    % 2. Calcul analytique de la dérivée dV/dt en injectant les équations du modèle
    % dV/dt = (1 - S*/S)*dS/dt + dE/dt + ((mu+alpha)/alpha)*dIdt
    dVdt = zeros(size(S));
    for k = 1:length(S)
        dS = Lambda - beta*S(k)*I(k) - mu*S(k);
        dE = beta*S(k)*I(k) - (mu + alpha)*E(k);
        dI = alpha*E(k) - (mu + gamma)*I(k);
        
        dVdt(k) = (1 - (S_star / S(k))) * dS + dE + ((mu + alpha) / alpha) * dI;
    end
end

% --- Équations du modèle SEIR ---
function dXdt = seir_eq(~, X, Lambda, beta, alpha, gamma, mu)
    S = X(1); E = X(2); I = X(3); R = X(4);
    dXdt = [Lambda - beta*S*I - mu*S; ...
            beta*S*I - (mu + alpha)*E; ...
            alpha*E - (mu + gamma)*I; ...
            gamma*I - mu*R];
end