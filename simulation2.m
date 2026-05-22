function main_SEIR()
    % --- Paramètres du modèle ---
    Lambda = 10;     % Taux de recrutement
    beta   = 0.002;  % Taux de transmission
    alpha  = 0.5;    % Taux de passage d'Exposé à Infecté
    gamma  = 0.3;    % Taux de guérison
    mu     = 0.1;    % Taux de mortalité naturelle

    % --- Conditions initiales ---
    S0 = 500; E0 = 10; I0 = 5; R0 = 0;
    X0 = [S0; E0; I0; R0]; 

    % --- Intervalle de temps ---
    tspan = [0 200];

    % --- Résolution Numérique ---
    [t, X] = ode45(@(t, X) SEIR_equations(t, X, Lambda, beta, alpha, gamma, mu), tspan, X0);

    % --- Tracé des Courbes (S, E, I, R) ---
    figure;
    plot(t, X(:,1), 'b-', 'LineWidth', 2); hold on;
    plot(t, X(:,2), 'r-', 'LineWidth', 2);
    plot(t, X(:,3), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2); % Jaune/Orange
    plot(t, X(:,4), 'm-', 'LineWidth', 2); % Violet
    hold off; grid on;
    
    xlabel('Temps (jours)', 'FontSize', 12);
    ylabel('Population', 'FontSize', 12);
    title('Simulation Numérique du Modèle SEIR via ode45', 'FontSize', 14);
    legend('Sains (S)', 'Exposés (E)', 'Infectés (I)', 'Guéris (R)', 'Location', 'best');
end

% --- Définition du champ de vecteurs ---
function dXdt = SEIR_equations(t, X, Lambda, beta, alpha, gamma, mu)
    S = X(1); E = X(2); I = X(3); R = X(4);
    dXdt = zeros(4, 1);
    
    dXdt(1) = Lambda - beta * S * I - mu * S;
    dXdt(2) = beta * S * I - (mu + alpha) * E;
    dXdt(3) = alpha * E - (mu + gamma) * I;
    dXdt(4) = gamma * I - mu * R;
end