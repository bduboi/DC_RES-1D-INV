function S = Sens_log(rho,h,AB)
%Calcule la matrice de sensitivité d'un modèle

% Calcul de la resitivité apparente à l'aide de l'opérateur direct
rhoa=dcfwdf((rho),(h),AB);

% Introduction d'une petite perturbation
delta=1e-3;

S    = zeros(length(rhoa),length(rho));
para = [rho; h];
for ii=1:length(rho) + length(h)
    
    % introduction des paramètres
    para_d=para;
para_d(ii)=(para(ii))*(1+delta); % Perturbation des paramètres
 % calcul de la réponse des données à des paramètres perturbés
rhoa_d = dcfwdf(para_d(1:length(rho)),para_d(length(rho)+1:end),AB);    
% Approximation de la matrice de sensibilité par la méthode des petites
% perturbations (schéma de dérivée du type [  f(x+h) - f(x) / h ]
S(:,ii) = (log(rhoa_d)-(log(rhoa)))./log(1+delta);
end