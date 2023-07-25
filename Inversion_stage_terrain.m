function [x, obj_fn, lambda_history] = Inversion_stage_terrain(x0,len_step,dobs,a,std)




%Paramètre de régularisation et de tolérance
lambda = 100; %paramètre de régularisation
lambda_droptol = 0.1; % De combien on tolère que lambda soit diminuer
lambda_scaling_factor = 0.1; % De combien on diminue lambda à chaque fois quand il y a convergence
lambda_min=1e-12; % valeur minimale que peut prendre lambda

tol=1e-4; % Valeur minimale tolérée pour le gradient
max_gn_iterations = 100; % itération max
c1 = 1e-4; % constante utilisée dans la recherche linéaire
obj_fn_goal=1e-8; % but souhaité pour la fonction objectif
goal_obj_diff = 1e-5; % but souhaité pour la différence de fonction objectif
step_size_min=1e-3; % valeur minimum du facteur appliqué à Delta m


%Define a starting model
%Definition des couches successives

I = eye(length(x0)); % Matrice identité
Weight_data = diag(1./std); % Matrice contenant les erreurs et appliquant un poids aux données
E = Weight_data'*Weight_data; % Quantité utilisée dans les équations, obtenue en calculant la norme de la fonction objectif modifée


% Fonction handles 
getData = @(x) dcfwdf(x(1:len_step),x(len_step+1:end),a); % getData(x) calcule directement les données sans avoir à rentrer a

residual = @(dobs, d) (log(dobs) - log(d));% Calcul la différence entre les données (log)

data_norm = @(r) 0.5 * norm(r)^2; % calcule la norme d'un vecteur

objective_function = @(r) data_norm(r); % même chose que data norm

%Transformation des paramètres pour éviter les problèmes de signe
log_x = log(x0); 

%Données générées par le modèle de départ

d = getData(x0);

%Résidu avec les données observées et calcul de la fonction objectif
Delta_d = residual(dobs,d);

obj_fn_current = data_norm(Delta_d); 

%%Initialisation des quantités pour garder la trace
obj_fn = zeros(max_gn_iterations, 1);
lambda_history = zeros(max_gn_iterations, 1);
norm_gc = Inf();
obj_fn_diff = Inf();

gn_it_counter = 1;

obj_fn(1) = obj_fn_current;

x=x0;
% Début de la boucle et définition de condition de terminaison
% indépendantes
figure('position', [0, 0, 500, 400])
while (norm_gc > tol) && ...
        (gn_it_counter < max_gn_iterations) && ...
        (obj_fn_current > obj_fn_goal) && ...
        (obj_fn_diff > goal_obj_diff)
    %Compute the Jacobian sensitivity matrix, computed without the log, no
    %chain rule

    J = Sens_log(x(1:len_step),x(len_step+1:end),a); %calcule la sensitivité
    gc = J'*E*Delta_d; % calcule le gradient
    norm_gc = norm(gc); % norme du gradient, utilisée dans les condition de terminaison

    %Calcule la direction du modèle à partir de l'équation normale
    Delta_x =  ((J'*E*J + lambda * I) \ gc);
    
    step_size = 1.0;
    log_x_proposed = log_x + step_size *Delta_x;
    d_proposed = getData(exp(log_x_proposed));
    
    obj_fn_proposed = objective_function(...
        residual(dobs, d_proposed));
    obj_fn_diff = obj_fn_current - obj_fn_proposed;

  % Recherche linéaire
    % derivée directionnelle pour la règle de Armijo dans la recherche
    % linéaire
    dir_deriv = -Delta_x' * gc;
    gn_it_counter = gn_it_counter + 1;
    %Début de la boucle de recherche linéaire. Diminue le facteur step_size
    %jusqu'à ce que obj_fn_diff redevienne supérieure au gradient * Delta x
    while (obj_fn_diff < -c1 * step_size * dir_deriv) % Doit être négatif pour assurer une descente
        step_size = step_size / 2.0;
        if step_size < step_size_min
            break%('Line search breakdown. Exit.');
        end
        log_x_proposed = log_x + step_size * Delta_x;
        d_proposed = getData(exp(log_x_proposed));
        obj_fn_proposed = objective_function(...
            residual(dobs, d_proposed));
        obj_fn_diff = obj_fn_current - obj_fn_proposed;
    end
    
    obj_rel_drop = obj_fn_diff / obj_fn_current;
    
    % En fin de boucle, mise à jour du paramètre x.
    x = exp(log_x_proposed);
    log_x = log_x_proposed;
    d_proposed = getData(x);% Calculer la nouvelle réponse du modèle (données)
    d = d_proposed;
    Delta_d = residual(dobs, d);
    obj_fn_current = obj_fn_proposed;
    
    % Faire décroitre le paramètre de régularisation dès qu'une certaine
    % convergence est atteinte (diminution relative de la fonction
    % objectif, fixé au début par obj_rel_drop-
    if obj_rel_drop < lambda_droptol && gn_it_counter > 1
        lambda = lambda * lambda_scaling_factor;
    end
    
    % Si lambda devient trop faible, il est mis à jour à une valeur
    % minimale décidée au début du code
    %
    if lambda < lambda_min
        lambda = lambda_min;
    end
    
    % Enregistrer l'évolution de lambda et de la fonction objectif.
    %
    lambda_history(gn_it_counter) = lambda;
    obj_fn(gn_it_counter) = obj_fn_current;
     % PLOT
    if mod(gn_it_counter, 1) == 0
        clf
        plot(a,dobs,'Color','blue','Marker','.','LineStyle',':','MarkerSize',4)
        hold on
        loglog(a,d,LineStyle="-",color='green');
        hold on
        errorbar(a,dobs,std,'vertical','Marker','|',"Color",'black','LineStyle','none', 'LineWidth', 0.5,'CapSize', 1)
        hold off
        hAx=gca;
        hAx.XScale='log';
        hAx.YScale='log';
        ylabel("Résistivité Apparente \rho_a (\Omega.m)")
        xlabel("AB/2 (m)")
        ylim([11,70])
        xlim([1,700])
        grid
        drawnow
        pause(0.2) 
    end
    
end
% Calcul du RMS. 
RMS = sqrt(sum(residual(dobs,d).^2 + std.^2));
legend("Resistivité Observée","Resistivité modélisée","Location",'northwest')
text(15,13," RMS = "+string(RMS),"Fontweight","bold")
end
