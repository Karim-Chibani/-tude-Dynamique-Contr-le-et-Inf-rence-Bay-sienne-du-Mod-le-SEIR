# Introduction Générale

La modélisation mathématique est devenue un outil indispensable pour comprendre, prédire et contrôler la propagation des maladies infectieuses au sein des populations. Parmi les modèles compartimentaux classiques, le modèle $SEIR$ se distingue par sa capacité à intégrer une phase d'incubation, essentielle pour l'étude de nombreuses pathologies contemporaines.

Pour le système $SEIR$ suivant :
$$
\left\{
\begin{aligned}
\frac{dS}{dt} &= \Lambda-\beta SI-\mu S \\
\frac{dE}{dt} &= \beta SI-(\mu+\alpha)E \\
\frac{dI}{dt} &= \alpha E-(\mu+\gamma)I \\
\frac{dR}{dt} &= \gamma I-\mu R
\end{aligned}
\right.
$$
Les conditions initiales représentent l’état de la population au temps initial $t = 0$. On écrit généralement :

$$S(0) = S_0, \quad E(0) = E_0, \quad I(0) = I_0, \quad R(0) = R_0$$
avec les contraintes de non-négativité :

$$S_0 \ge 0, \quad E_0 \ge 0, \quad I_0 \ge 0, \quad R_0 \ge 0$$

et la condition sur la population totale initiale :

$$N_0 = S_0 + E_0 + I_0 + R_0 > 0$$

où :
* **$S_0$** : Nombre initial d’individus susceptibles.
* **$E_0$** : Nombre initial d’individus exposés.
* **$I_0$** : Nombre initial d’individus infectés.
* **$R_0$** : Nombre initial d’individus rétablis.
où $S$, $E$, $I$, et $R$ représentent respectivement les fractions de la population Susceptible, Exposée (en incubation), Infectée et Retirée (guérie ou immunisée).
l'analyse mathématique classique permet de caractériser la dynamique locale et globale de l'épidémie à travers le taux de reproduction de base $R_0$. Cependant, l'application de ces modèles théoriques à des scénarios réels soulève deux défis majeurs : l'élaboration de stratégies d'intervention optimales et l'identification précise des paramètres biologiques à partir de données de terrain souvent incomplètes ou bruitées.

Ce projet s'articule autour de ces défis à travers une double approche.

J'ai intégré vos choix de simulation ($\beta = 0.002$ et $\beta = 0.008$) et j'ai enrichi le texte avec les concepts mathématiques rigoureux (le point d'équilibre sans maladie $DFE$ et le point d'équilibre endémique $EE$) pour donner une dimension hautement scientifique et professionnelle à votre travail.

nous développons un cadre de contrôle optimal appliqué à ce système d'équations différentielles ordinaires afin de planifier des politiques sanitaires (telles que la vaccination ou le confinement) capables de minimiser l'impact épidémique à moindre coût.

nous abordons le problème inverse de l'estimation paramétrique. Les données collectées lors des crises sanitaires étant intrinsèquement entachées d'erreurs de mesure et de fluctuations aléatoires, nous mettons en œuvre une approche statistique robuste fondée sur l'inférence bayésienne. 

En utilisant des techniques d'échantillonnage de Monte Carlo par Chaînes de Markov (MCMC), notamment l'algorithme de Metropolis-Hastings, ce travail vise à filtrer le bruit numérique pour reconstruire les distributions de probabilité *a posteriori* des paramètres clés, tels que le taux de transmission $\beta$, garantissant ainsi des prévisions épidémiologiques hautement fiables et scientifiquement rigoureuses.

### 1. Le cadre fonctionnel

Comme il s’agit d’un système d’équations différentielles ordinaires (EDO), le cadre fonctionnel naturel est un espace de fonctions continues ou continûment dérivables.

On pose le vecteur d'état :

$$X(t) = \big(S(t), E(t), I(t), R(t)\big)$$

Alors, l'application $X$ est définie sur un horizon temporel $[0, T]$ (avec $T > 0$) par :

$$X : [0, T] \to \mathbb{R}^4$$

Le cadre fonctionnel classique pour garantir l'existence et l'unicité des solutions (via le théorème de Cauchy-Lipschitz) est l'espace des fonctions de classe $C^1$ :

$$X \in C^1([0, T], \mathbb{R}^4)$$

C’est-à-dire que :
* Les composantes $S, E, I, R$ sont continues sur $[0, T]$.
* Leurs dérivées temporelles sont également continues sur $[0, T]$.
  ### 2. Existence locale et unicité

Pour établir l'existence locale et l'unicité des solutions, on applique le **théorème de Cauchy-Lipschitz** (parfois appelé théorème de Picard-Lindelöf). La condition fondamentale à vérifier est la régularité du champ de vecteurs $F$, plus précisément son caractère localement lipschitzien.

---

### 3. Régularité de $F$

Les composantes de l'application $F = (F_1, F_2, F_3, F_4)$ sont définies par :

$$\begin{aligned}
F_1(S,E,I,R) &= \Lambda - \beta SI - \mu S \\
F_2(S,E,I,R) &= \beta SI - (\mu + \alpha)E \\
F_3(S,E,I,R) &= \alpha E - (\mu + \gamma)I \\
F_4(S,E,I,R) &= \gamma I - \mu R
\end{aligned}$$

Chacune de ces composantes est une fonction polynomiale par rapport aux variables $S, E, I, R$. Par conséquent, le champ de vecteurs est de classe $C^\infty$ sur $\mathbb{R}^4$ :

$$F \in C^\infty(\mathbb{R}^4, \mathbb{R}^4)$$

Puisque toute fonction de classe $C^1$ (et *a fortiori* de classe $C^\infty$) est localement lipschitzienne sur son domaine de définition, on en conclut directement que **$F$ est localement lipschitzienne**.

---

### 4. Conclusion locale

Par application directe du théorème de Cauchy-Lipschitz, pour toute condition initiale $X(0) = X_0 \in \mathbb{R}^4$, il existe un temps d'existence maximal $T_{\max} > 0$ et une unique solution maximale :

$$X \in C^1([0, T_{\max}), \mathbb{R}^4)$$

qui satisfait le problème de Cauchy associé.

---

### 5. Positivité et invariance des solutions

Pour que le modèle soit biologiquement cohérent, les variables d'état représentant des effectifs de population doivent rester non négatives pour tout temps $t \ge 0$. On montre que le domaine de réalité biologique $\mathbb{R}_+^4$ est **positivement invariant** en analysant les signes des dérivées sur les bords de ce domaine :

* **Pour $S(t)$ :** Si $S = 0$ (avec $E, I, R \ge 0$), alors $\frac{dS}{dt} = \Lambda \ge 0$.
* **Pour $E(t)$ :** Si $E = 0$ (avec $S, I, R \ge 0$), alors $\frac{dE}{dt} = \beta SI \ge 0$.
* **Pour $I(t)$ :** Si $I = 0$ (avec $S, E, R \ge 0$), alors $\frac{dI}{dt} = \alpha E \ge 0$.
* **Pour $R(t)$ :** Si $R = 0$ (avec $S, E, I \ge 0$), alors $\frac{dR}{dt} = \gamma I \ge 0$.

Puisque le flux du champ de vecteurs $F$ pointe vers l'intérieur (ou reste tangent) aux frontières de l'hyperoctant positif, aucune trajectoire issue de $\mathbb{R}_+^4$ ne peut franchir ces frontières.

> **Conclusion :** L'ensemble $\mathbb{R}_+^4$ est un sous-ensemble positivement invariant pour le système. Si $X_0 \in \mathbb{R}_+^4$, alors $X(t) \in \mathbb{R}_+^4$ pour tout $t \in [0, T_{\max})$.

---

### 6. Existence globale

Pour étendre la solution locale en une solution globale (c'est-à-dire prouver que $T_{\max} = +\infty$), il suffit d'exclure tout phénomène d'explosion en temps fini. On introduit pour cela la population totale $N(t)$ définie par :

$$N(t) = S(t) + E(t) + I(t) + R(t)$$

---

### 7. Équation sur la population totale

En sommant les quatre équations différentielles du système SEIR, les termes d'interaction de masse se simplifient deux à deux ($\beta SI$, $\alpha E$ et $\gamma I$). On obtient l'EDO linéaire suivante sur $N(t)$ :

$$\frac{dN}{dt} = \Lambda - \mu S - \mu E - \mu I - \mu R = \Lambda - \mu N$$

---

### 8. Bornes a priori sur la solution

La résolution de cette équation différentielle linéaire du premier ordre donne explicitement :

$$N(t) = N_0 e^{-\mu t} + \frac{\Lambda}{\mu} \left(1 - e^{-\mu t}\right)$$

En étudiant le comportement asymptotique ou en appliquant un principe de comparaison, on en déduit que la population totale est uniformément bornée supérieurement :

$$0 \le N(t) \le \max\left(N_0, \frac{\Lambda}{\mu}\right)$$

Comme toutes les variables d'état sont non négatives ($S, E, I, R \ge 0$), chacune d'elles est trivialement majorée par la population totale $N(t)$. Par conséquent, les composantes du vecteur d'état $X(t)$ restent confinées dans un sous-ensemble compact de $\mathbb{R}_+^4$.

---

### 9. Conclusion globale

En résumé :
1. Une unique solution locale existe sur l'intervalle $[0, T_{\max})$.
2. Le champ de vecteurs $F$ est localement lipschitzien.
3. Les trajectoires de la solution restent bornées dans un compact de $\mathbb{R}_+^4$.

D'après le **théorème de prolongement des solutions** (alternative de l'explosion), la solution ne peut pas exploser en temps fini, ce qui implique nécessairement que le temps d'existence maximal est infini ($T_{\max} = +\infty$).

La solution du problème de Cauchy est donc globale et s'écrit :

$$X \in C^1\left([0, +\infty), \mathbb{R}_+^4\right)$$
